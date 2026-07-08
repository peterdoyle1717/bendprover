/* embcheck_mp.c -- certified embeddedness checker for bend records.
 *
 * Mechanical port of embcheck.py/mpball.py to C-MPFR (the python pair
 * stays as the differential-test oracle). Consumes FORMAT.md records
 * (halfturn units, integer tokens exact, benderr per record) and
 * certifies, per net:
 *
 *   EDGE   shared-edge faces cross iff |bend| = 1 halfturn; certify
 *          |b| + benderr < 1 (integer tokens +-1 reject immediately).
 *   LINK   faces sharing exactly a vertex cross near it iff the vertex
 *          link polygon self-intersects on S^2; runs of exactly-flat
 *          corners merge (coplanar arcs share a great circle); lune
 *          (2 runs) is covered by EDGE; otherwise non-adjacent merged
 *          arcs are certified disjoint (per-sub-arc hemisphere fast
 *          path, else exclusion of the circle intersection points).
 *   PAIR   vertex-disjoint simplex pairs (V+E+F): double-precision
 *          bounding-ball cull with a conservative margin, then
 *          candidate-axis separation in ball arithmetic.
 *
 * Ball arithmetic: 3-vector balls (mpfr centers at BPREC bits, one
 * double radius), radius propagation with center norms over-estimated
 * in double and an ulp term per operation. Enclosures contain the
 * exact realization's vertices by construction of the development.
 *
 * LEMMA DEBT (as in embcheck.py): the shared-vertex reduction and the
 * per-op inflation constants await written proof / independent review.
 *
 * Build: cc -O2 -o embcheck_mp embcheck_mp.c -lmpfr -lgmp
 * Usage: embcheck_mp RECORDS.bends...   (TSV verdicts on stdout)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpfr.h>

#define BPREC 192
#define MAXV 128
#define MAXF (2 * MAXV + 4)
#define MAXE (3 * MAXV + 6)
#define MAXDEG 16
#define MAXS (MAXV + MAXE + MAXF)
#define ULPC ldexp(1.0, 20 - BPREC)
#define COS_ARC 0.5           /* sub-arcs have length pi/3 */
#define CULL_MARGIN 0.05

static int NV, NF, NE;
static int FACES[MAXF][3];
static int EDGE_A[MAXE], EDGE_B[MAXE];
static int EDGE_IDX[MAXV + 1][MAXV + 1];
static int APEX[MAXV + 1][MAXV + 1];      /* directed (a,b) -> third vertex */

/* ---------------- scalar and vector balls ---------------- */

typedef struct { mpfr_t c; double r; } SB;
typedef struct { mpfr_t c[3]; double r; } VB;

static void sb_init(SB *s) { mpfr_init2(s->c, BPREC); s->r = 0; }
static void vb_init(VB *v) { for (int k = 0; k < 3; k++) mpfr_init2(v->c[k], BPREC); v->r = 0; }

static double dnorm(const VB *v)          /* over-estimate of |center| */
{
    double s = 0;
    for (int k = 0; k < 3; k++) { double x = mpfr_get_d(v->c[k], MPFR_RNDA); s += x * x; }
    return sqrt(s) * (1 + 1e-13) + 1e-300;
}
static double dabs(const SB *s) { return fabs(mpfr_get_d(s->c, MPFR_RNDA)) * (1 + 1e-13) + 1e-300; }

static void vb_sub(VB *o, const VB *a, const VB *b)
{
    for (int k = 0; k < 3; k++) mpfr_sub(o->c[k], a->c[k], b->c[k], MPFR_RNDN);
    o->r = a->r + b->r; o->r += ULPC * (dnorm(o) + 1);
}
static void vb_add(VB *o, const VB *a, const VB *b)
{
    for (int k = 0; k < 3; k++) mpfr_add(o->c[k], a->c[k], b->c[k], MPFR_RNDN);
    o->r = a->r + b->r; o->r += ULPC * (dnorm(o) + 1);
}
static void vb_scale_sb(VB *o, const SB *s, const VB *v)   /* o = s*v, o != v ok aliasing avoided by caller */
{
    double nv = dnorm(v), ns = dabs(s);
    for (int k = 0; k < 3; k++) mpfr_mul(o->c[k], v->c[k], s->c, MPFR_RNDN);
    o->r = ns * v->r + nv * s->r + s->r * v->r;
    o->r = o->r * (1 + 1e-13) + ULPC * (dnorm(o) + 1);
}
static void vb_scale_d(VB *o, double d, const VB *v)       /* exact double factor */
{
    for (int k = 0; k < 3; k++) mpfr_mul_d(o->c[k], v->c[k], d, MPFR_RNDN);
    o->r = fabs(d) * v->r * (1 + 1e-13) + ULPC * (dnorm(o) + 1);
}
static void sb_dot(SB *o, const VB *a, const VB *b)
{
    static mpfr_t t; static int init = 0;
    if (!init) { mpfr_init2(t, BPREC); init = 1; }
    mpfr_mul(o->c, a->c[0], b->c[0], MPFR_RNDN);
    mpfr_mul(t, a->c[1], b->c[1], MPFR_RNDN); mpfr_add(o->c, o->c, t, MPFR_RNDN);
    mpfr_mul(t, a->c[2], b->c[2], MPFR_RNDN); mpfr_add(o->c, o->c, t, MPFR_RNDN);
    double na = dnorm(a), nb = dnorm(b);
    o->r = (na * b->r + nb * a->r + a->r * b->r) * (1 + 1e-13) + ULPC * (dabs(o) + 1);
}
static void vb_cross(VB *o, const VB *a, const VB *b)      /* o distinct from a,b */
{
    static mpfr_t t; static int init = 0;
    if (!init) { mpfr_init2(t, BPREC); init = 1; }
    for (int k = 0; k < 3; k++) {
        int i = (k + 1) % 3, j = (k + 2) % 3;
        mpfr_mul(o->c[k], a->c[i], b->c[j], MPFR_RNDN);
        mpfr_mul(t, a->c[j], b->c[i], MPFR_RNDN);
        mpfr_sub(o->c[k], o->c[k], t, MPFR_RNDN);
    }
    double na = dnorm(a), nb = dnorm(b);
    o->r = (na * b->r + nb * a->r + a->r * b->r) * (1 + 1e-13) + ULPC * (dnorm(o) + 1);
}
static int vb_unit(VB *o, const VB *v)    /* o = v/|v|; -1 if not certified > 0 */
{
    static SB n2, n; static int init = 0;
    if (!init) { sb_init(&n2); sb_init(&n); init = 1; }
    sb_dot(&n2, v, v);
    double lo = mpfr_get_d(n2.c, MPFR_RNDD) - n2.r;
    if (!(lo > 0)) return -1;
    mpfr_sqrt(n.c, n2.c, MPFR_RNDN);
    n.r = n2.r / (2 * sqrt(lo)) * (1 + 1e-13) + ULPC;
    double nlo = mpfr_get_d(n.c, MPFR_RNDD) - n.r;
    if (!(nlo > 0)) return -1;
    double nv = dnorm(v);
    for (int k = 0; k < 3; k++) mpfr_div(o->c[k], v->c[k], n.c, MPFR_RNDN);
    o->r = (v->r + nv * n.r / nlo) / nlo * (1 + 1e-13) + ULPC * 2;
    return 0;
}
static double sb_lo(const SB *s) { return mpfr_get_d(s->c, MPFR_RNDD) - s->r; }
static double sb_hi(const SB *s) { return mpfr_get_d(s->c, MPFR_RNDU) + s->r; }

/* ---------------- record state ---------------- */

static char NETNAME[256];
static double BENDERR;                    /* halfturns */
static char *BVAL[MAXE];                  /* value token per edge, malloc'd */
static VB POS[MAXV + 1];
static int POS_INIT = 0;

/* bend angle ball in radians for edge e */
static void bend_ball(SB *th, int e)
{
    static mpfr_t pi; static int init = 0;
    if (!init) { mpfr_init2(pi, BPREC); mpfr_const_pi(pi, MPFR_RNDN); init = 1; }
    const char *s = BVAL[e];
    if (!strcmp(s, "0")) { mpfr_set_ui(th->c, 0, MPFR_RNDN); th->r = 0; return; }
    if (!strcmp(s, "1"))  { mpfr_set(th->c, pi, MPFR_RNDN); th->r = ULPC; return; }
    if (!strcmp(s, "-1")) { mpfr_neg(th->c, pi, MPFR_RNDN); th->r = ULPC; return; }
    mpfr_set_str(th->c, s, 10, MPFR_RNDN);
    mpfr_mul(th->c, th->c, pi, MPFR_RNDN);
    th->r = BENDERR * 3.14159265358979324 * (1 + 1e-13) + ULPC * 4;
}

/* ---------------- topology from the faces line ---------------- */

static int build_topology(const char *nc)
{
    NV = 0; NF = 0; NE = 0;
    const char *p = nc;
    while (*p) {
        int a, b, c;
        if (sscanf(p, "%d,%d,%d", &a, &b, &c) != 3) return -1;
        if (NF >= MAXF || a < 1 || b < 1 || c < 1 || a > MAXV || b > MAXV || c > MAXV) return -1;
        FACES[NF][0] = a; FACES[NF][1] = b; FACES[NF][2] = c; NF++;
        if (a > NV) NV = a;
        if (b > NV) NV = b;
        if (c > NV) NV = c;
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }
    for (int i = 0; i <= NV; i++)
        for (int j = 0; j <= NV; j++) { EDGE_IDX[i][j] = -1; APEX[i][j] = 0; }
    for (int f = 0; f < NF; f++) {
        int v[3] = { FACES[f][0], FACES[f][1], FACES[f][2] };
        for (int k = 0; k < 3; k++) {
            int a = v[k], b = v[(k + 1) % 3], c = v[(k + 2) % 3];
            APEX[a][b] = c;
            int lo = a < b ? a : b, hi = a < b ? b : a;
            if (EDGE_IDX[lo][hi] < 0) {
                if (NE >= MAXE) return -1;
                EDGE_IDX[lo][hi] = NE; EDGE_A[NE] = lo; EDGE_B[NE] = hi; NE++;
            }
        }
    }
    return 0;
}

/* ---------------- ball development (place_third offset form) ---------------- */

static int develop(void)
{
    static int placedf[MAXF], queue[MAXF];
    static VB M, E, W, N, T1, T2; static SB CT, ST, TH, S32;
    static int init = 0;
    if (!init) {
        vb_init(&M); vb_init(&E); vb_init(&W); vb_init(&N); vb_init(&T1); vb_init(&T2);
        sb_init(&CT); sb_init(&ST); sb_init(&TH); sb_init(&S32);
        mpfr_sqrt_ui(S32.c, 3, MPFR_RNDN);
        mpfr_div_ui(S32.c, S32.c, 2, MPFR_RNDN);
        S32.r = ULPC;                       /* sqrt(3)/2 at BPREC bits */
        init = 1;
    }
    if (!POS_INIT) { for (int v = 0; v <= MAXV; v++) vb_init(&POS[v]); POS_INIT = 1; }
    int b0 = FACES[0][1], c0 = FACES[0][2];
    for (int v = 0; v <= NV; v++) { POS[v].r = 0; for (int k = 0; k < 3; k++) mpfr_set_ui(POS[v].c[k], 0, MPFR_RNDN); }
    mpfr_set_ui(POS[b0].c[0], 1, MPFR_RNDN);
    mpfr_set_d(POS[c0].c[0], 0.5, MPFR_RNDN);
    mpfr_sqrt_ui(POS[c0].c[1], 3, MPFR_RNDN);
    mpfr_div_ui(POS[c0].c[1], POS[c0].c[1], 2, MPFR_RNDN);
    POS[c0].r = ULPC;
    for (int f = 0; f < NF; f++) placedf[f] = 0;
    placedf[0] = 1;
    int qh = 0, qt = 0; queue[qt++] = 0;
    /* face id lookup for directed edge via APEX + a face map */
    static int FID[MAXV + 1][MAXV + 1];
    for (int i = 0; i <= NV; i++) for (int j = 0; j <= NV; j++) FID[i][j] = -1;
    for (int f = 0; f < NF; f++) {
        int *v = FACES[f];
        FID[v[0]][v[1]] = f; FID[v[1]][v[2]] = f; FID[v[2]][v[0]] = f;
    }
    while (qh < qt) {
        int f = queue[qh++];
        int v[3] = { FACES[f][0], FACES[f][1], FACES[f][2] };
        for (int k = 0; k < 3; k++) {
            int x = v[k], y = v[(k + 1) % 3], p = v[(k + 2) % 3];
            int g = FID[y][x];
            if (g < 0 || placedf[g]) continue;
            int d = APEX[y][x];
            /* M = (X+Y)/2; w = unit(P-M); n = unit(w x (Y-X));
               D = M + (sqrt3/2)(-cos th * w + sin th * n) */
            vb_add(&T1, &POS[x], &POS[y]); vb_scale_d(&M, 0.5, &T1);
            vb_sub(&E, &POS[y], &POS[x]);
            vb_sub(&T1, &POS[p], &M);
            if (vb_unit(&W, &T1) < 0) {
                if (getenv("EMB_DEBUG"))
                    fprintf(stderr, "unit-W fail placing %d from face %d: T1 r=%.3g norm=%.3g posr=%.3g/%.3g/%.3g\n",
                            d, f, T1.r, dnorm(&T1), POS[x].r, POS[y].r, POS[p].r);
                return -1;
            }
            vb_cross(&T1, &W, &E);
            if (vb_unit(&N, &T1) < 0) {
                if (getenv("EMB_DEBUG"))
                    fprintf(stderr, "unit-N fail placing %d from face %d: T1 r=%.3g norm=%.3g W.r=%.3g E.r=%.3g\n",
                            d, f, T1.r, dnorm(&T1), W.r, E.r);
                return -1;
            }
            int lo = x < y ? x : y, hi = x < y ? y : x;
            bend_ball(&TH, EDGE_IDX[lo][hi]);
            mpfr_cos(CT.c, TH.c, MPFR_RNDN); CT.r = TH.r + ULPC;
            mpfr_neg(CT.c, CT.c, MPFR_RNDN);
            mpfr_sin(ST.c, TH.c, MPFR_RNDN); ST.r = TH.r + ULPC;
            vb_scale_sb(&T1, &CT, &W);
            vb_scale_sb(&T2, &ST, &N);
            vb_add(&T1, &T1, &T2);
            vb_scale_sb(&T2, &S32, &T1);
            vb_add(&POS[d], &M, &T2);
            if (getenv("EMB_TRACE"))
                fprintf(stderr, "place %d r=%.6g (x=%d rx=%.3g y=%d ry=%.3g p=%d rp=%.3g W.r=%.3g N.r=%.3g)\n",
                        d, POS[d].r, x, POS[x].r, y, POS[y].r, p, POS[p].r, W.r, N.r);
            placedf[g] = 1; queue[qt++] = g;
        }
    }
    return 0;
}

/* ---------------- link simplicity with merged flat runs ---------------- */

static int link_check(char *why, size_t wl)
{
    static VB DIR[MAXDEG], T; static int init = 0;
    static VB NP, NQ, XR, X;
    if (!init) {
        for (int i = 0; i < MAXDEG; i++) vb_init(&DIR[i]);
        vb_init(&T); vb_init(&NP); vb_init(&NQ); vb_init(&XR); vb_init(&X);
        init = 1;
    }
    static SB S, S1, S2; static int sinit = 0;
    if (!sinit) { sb_init(&S); sb_init(&S1); sb_init(&S2); sinit = 1; }
    for (int v = 1; v <= NV; v++) {
        /* cyclic star */
        int cyc[MAXDEG], d = 0, w0 = -1;
        for (int w = 1; w <= NV && w0 < 0; w++) if (w != v && APEX[v][w]) w0 = w;
        if (w0 < 0) continue;
        int w = w0;
        do {
            if (d >= MAXDEG) return -1;
            cyc[d++] = w;
            w = APEX[v][w];
        } while (w != w0);
        int flatc[MAXDEG], anyfold = 0;
        for (int i = 0; i < d; i++) {
            int a = v < cyc[i] ? v : cyc[i], b = v < cyc[i] ? cyc[i] : v;
            flatc[i] = !strcmp(BVAL[EDGE_IDX[a][b]], "0");
            if (!flatc[i]) anyfold = 1;
        }
        if (!anyfold) continue;                    /* fully flat vertex */
        for (int i = 0; i < d; i++) {
            vb_sub(&T, &POS[cyc[i]], &POS[v]);
            if (vb_unit(&DIR[i], &T) < 0) { snprintf(why, wl, "degenerate direction at %d", v); return 1; }
        }
        /* merged runs: split at non-flat corners */
        int start = -1;
        for (int i = 0; i < d; i++) if (!flatc[i]) { start = i; break; }
        int runs[MAXDEG][MAXDEG + 1], rlen[MAXDEG], m = 0, cur = 0;
        runs[0][0] = start; cur = 1;
        for (int k = 1; k <= d; k++) {
            int i = (start + k) % d;
            runs[m][cur++] = i;
            if (!flatc[i]) { rlen[m++] = cur; runs[m][0] = i; cur = 1; }
        }
        if (m == 1) { snprintf(why, wl, "single-fold link at vertex %d", v); return 1; }
        if (m == 2) continue;                      /* lune; covered by EDGE */
        for (int i = 0; i < m; i++) {
            for (int j = i + 1; j < m; j++) {
                int dj = j - i;
                if (dj == 1 || dj == m - 1) continue;   /* adjacent merged arcs */
                /* fast path: all of Q one side of circle(P) (per-sub-arc lemma) */
                vb_cross(&NP, &DIR[runs[i][0]], &DIR[runs[i][1]]);
                vb_cross(&NQ, &DIR[runs[j][0]], &DIR[runs[j][1]]);
                int okpair = 0;
                for (int side = 0; side < 2 && !okpair; side++) {
                    const VB *nrm = side ? &NQ : &NP;
                    int ri = side ? i : j;
                    int allpos = 1, allneg = 1;
                    for (int t = 0; t < rlen[ri]; t++) {
                        sb_dot(&S, nrm, &DIR[runs[ri][t]]);
                        if (!(sb_lo(&S) > 0)) allpos = 0;
                        if (!(sb_hi(&S) < 0)) allneg = 0;
                    }
                    if (allpos || allneg) okpair = 1;
                }
                if (okpair) continue;
                /* exclusion of the circle intersection points +-x */
                vb_cross(&XR, &NP, &NQ);
                sb_dot(&S, &XR, &XR);
                if (!(sb_lo(&S) > 0)) {
                    snprintf(why, wl, "link arcs cocircular at vertex %d", v); return 1;
                }
                if (vb_unit(&X, &XR) < 0) { snprintf(why, wl, "link x-unit at vertex %d", v); return 1; }
                for (int sgn = 0; sgn < 2; sgn++) {
                    if (sgn) for (int k2 = 0; k2 < 3; k2++) mpfr_neg(X.c[k2], X.c[k2], MPFR_RNDN);
                    int outP = 1, outQ = 1;
                    for (int t = 0; t + 1 < rlen[i]; t++) {
                        sb_dot(&S1, &X, &DIR[runs[i][t]]);
                        sb_dot(&S2, &X, &DIR[runs[i][t + 1]]);
                        if (!(sb_hi(&S1) < COS_ARC || sb_hi(&S2) < COS_ARC)) { outP = 0; break; }
                    }
                    for (int t = 0; t + 1 < rlen[j]; t++) {
                        sb_dot(&S1, &X, &DIR[runs[j][t]]);
                        sb_dot(&S2, &X, &DIR[runs[j][t + 1]]);
                        if (!(sb_hi(&S1) < COS_ARC || sb_hi(&S2) < COS_ARC)) { outQ = 0; break; }
                    }
                    if (!outP && !outQ) {
                        snprintf(why, wl, "link not certified simple at vertex %d", v);
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}

/* ---------------- pair separation with double cull ---------------- */

static int SIMP[MAXS][3], SN[MAXS], NS;
static double SC[MAXS][3], SR[MAXS];

static void build_simplices(void)
{
    NS = 0;
    for (int v = 1; v <= NV; v++) { SIMP[NS][0] = v; SN[NS] = 1; NS++; }
    for (int e = 0; e < NE; e++) { SIMP[NS][0] = EDGE_A[e]; SIMP[NS][1] = EDGE_B[e]; SN[NS] = 2; NS++; }
    for (int f = 0; f < NF; f++) {
        SIMP[NS][0] = FACES[f][0]; SIMP[NS][1] = FACES[f][1]; SIMP[NS][2] = FACES[f][2];
        SN[NS] = 3; NS++;
    }
    for (int s = 0; s < NS; s++) {
        double c[3] = {0, 0, 0}, rmax = 0;
        for (int k = 0; k < SN[s]; k++)
            for (int j = 0; j < 3; j++)
                c[j] += mpfr_get_d(POS[SIMP[s][k]].c[j], MPFR_RNDN) / SN[s];
        for (int k = 0; k < SN[s]; k++) {
            double dd = 0;
            for (int j = 0; j < 3; j++) {
                double x = mpfr_get_d(POS[SIMP[s][k]].c[j], MPFR_RNDN) - c[j];
                dd += x * x;
            }
            double rr = sqrt(dd) + POS[SIMP[s][k]].r + 1e-12;
            if (rr > rmax) rmax = rr;
        }
        for (int j = 0; j < 3; j++) SC[s][j] = c[j];
        SR[s] = rmax;
    }
}

static int share_vertex(int a, int b)
{
    for (int i = 0; i < SN[a]; i++)
        for (int j = 0; j < SN[b]; j++)
            if (SIMP[a][i] == SIMP[b][j]) return 1;
    return 0;
}

static int pair_check(char *why, size_t wl)
{
    static VB AX; static SB S; static int init = 0;
    if (!init) { vb_init(&AX); sb_init(&S); init = 1; }
    build_simplices();
    for (int a = 0; a < NS; a++) {
        for (int b = a + 1; b < NS; b++) {
            if (share_vertex(a, b)) continue;
            double d2 = 0;
            for (int j = 0; j < 3; j++) { double x = SC[a][j] - SC[b][j]; d2 += x * x; }
            if (sqrt(d2) > SR[a] + SR[b] + CULL_MARGIN) continue;   /* culled, sound */
            /* candidate axes: centroid diff (exact double), vertex diffs */
            int ok = 0;
            for (int cand = 0; cand <= SN[a] * SN[b] && !ok; cand++) {
                double ax[3];
                if (cand == 0) { for (int j = 0; j < 3; j++) ax[j] = SC[a][j] - SC[b][j]; }
                else {
                    int i = (cand - 1) / SN[b], j2 = (cand - 1) % SN[b];
                    for (int j = 0; j < 3; j++)
                        ax[j] = mpfr_get_d(POS[SIMP[a][i]].c[j], MPFR_RNDN)
                              - mpfr_get_d(POS[SIMP[b][j2]].c[j], MPFR_RNDN);
                }
                double n = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
                if (n < 1e-12) continue;
                for (int j = 0; j < 3; j++) mpfr_set_d(AX.c[j], ax[j], MPFR_RNDN);
                AX.r = 0;                                       /* chosen exact axis */
                double loA = 1e300, hiA = -1e300, loB = 1e300, hiB = -1e300;
                for (int i = 0; i < SN[a]; i++) {
                    sb_dot(&S, &AX, &POS[SIMP[a][i]]);
                    if (sb_lo(&S) < loA) loA = sb_lo(&S);
                    if (sb_hi(&S) > hiA) hiA = sb_hi(&S);
                }
                for (int i = 0; i < SN[b]; i++) {
                    sb_dot(&S, &AX, &POS[SIMP[b][i]]);
                    if (sb_lo(&S) < loB) loB = sb_lo(&S);
                    if (sb_hi(&S) > hiB) hiB = sb_hi(&S);
                }
                if (loA > hiB || loB > hiA) ok = 1;
            }
            if (!ok) {
                snprintf(why, wl, "separation not certified (simplices %d,%d)", a, b);
                return 1;
            }
        }
    }
    return 0;
}

/* ---------------- record loop ---------------- */

/* returns 0 PASS, 1 FAIL, 2 PANCAKE (exactly flat, certified by exactness) */
static int check_record(const char *nc, char *why, size_t wl)
{
    if (build_topology(nc) < 0) { snprintf(why, wl, "bad faces"); return 1; }
    int all_int = 1;
    for (int e = 0; e < NE; e++) {
        if (!BVAL[e]) { snprintf(why, wl, "missing bend"); return 1; }
        const char *s = BVAL[e];
        if (strcmp(s, "0") && strcmp(s, "1") && strcmp(s, "-1")) all_int = 0;
    }
    if (all_int) {
        snprintf(why, wl, "all bends integer halfturns (flat pancake, not embedded)");
        return 2;
    }
    for (int e = 0; e < NE; e++) {
        const char *s = BVAL[e];
        if (!strcmp(s, "1") || !strcmp(s, "-1") ||
            (strcmp(s, "0") && fabs(atof(s)) + BENDERR >= 1.0 - 1e-14)) {
            snprintf(why, wl, "edge (%d,%d) bend at gem/2", EDGE_A[e], EDGE_B[e]);
            return 1;
        }
    }
    if (develop() < 0) { snprintf(why, wl, "development failed"); return 1; }
    int rc = link_check(why, wl);
    if (rc) return rc;
    return pair_check(why, wl);
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s RECORDS...\n", argv[0]); return 2; }
    printf("net\tv\tembedded\twhy\n");
    char line[262144], nc[262144], why[256];
    int npass = 0, nfail = 0;
    for (int ai = 1; ai < argc; ai++) {
        FILE *fp = fopen(argv[ai], "r");
        if (!fp) { fprintf(stderr, "cannot open %s\n", argv[ai]); return 2; }
        NETNAME[0] = '\0'; nc[0] = '\0'; BENDERR = 0;
        for (int e = 0; e < MAXE; e++) { free(BVAL[e]); BVAL[e] = NULL; }
        while (fgets(line, sizeof line, fp)) {
            if (!strncmp(line, "net ", 4)) {
                sscanf(line + 4, "%255s", NETNAME);
                nc[0] = '\0'; BENDERR = 0;
                for (int e = 0; e < MAXE; e++) { free(BVAL[e]); BVAL[e] = NULL; }
            } else if (!strncmp(line, "benderr ", 8)) {
                BENDERR = atof(line + 8);
            } else if (!strncmp(line, "faces ", 6)) {
                sscanf(line + 6, "%262000s", nc);
                if (build_topology(nc) < 0) nc[0] = '\0';
            } else if (!strncmp(line, "b ", 2)) {
                int a, b; char val[128];
                if (sscanf(line + 2, "%d %d %127s", &a, &b, val) == 3 &&
                    a >= 1 && b >= 1 && a <= MAXV && b <= MAXV && EDGE_IDX[a][b] >= 0)
                    BVAL[EDGE_IDX[a][b]] = strdup(val);
            } else if (!strncmp(line, "end", 3) && nc[0]) {
                why[0] = '\0';
                int rc = check_record(nc, why, sizeof why);
                printf("%s\t%d\t%s\t%s\n", NETNAME, NV,
                       rc == 0 ? "PASS" : rc == 2 ? "PANCAKE" : "FAIL", why);
                if (rc == 1) nfail++; else npass++;
                nc[0] = '\0';
            }
        }
        fclose(fp);
    }
    fprintf(stderr, "done: %d embedded-or-pancake, %d FAIL\n", npass, nfail);
    return 0;
}

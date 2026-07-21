/* euclid_lm_mp.c -- the LM solver of euclid_clean.c, converted to MPFR.
 *
 * Mechanical port: same pipeline (netcode -> topology -> wish init ->
 * LM-with-dent-gate at alpha = pi/3 -> realize -> OBJ), same LM logic
 * and constants. Differences from the DP original:
 *   - all solver arithmetic in MPFR at --prec bits (default 128);
 *   - the LM normal solve uses an MPFR Cholesky of (J^T J + lambda D)
 *     instead of LAPACK dgesv (the matrix is SPD);
 *   - the normal matrix is assembled by the vertex loop (each Jacobian
 *     column is supported on its edge's two flowers), same values as
 *     the dense dot products;
 *   - wish init solves its small V x V system in plain doubles (it is
 *     a starting guess; no LAPACK dependency);
 *   - OBJ output carries "# bend e a b value" comment lines and a
 *     "# resid" line, coordinates at %.40Rf.
 *
 * Build:  cc -O2 -o euclid_lm_mp euclid_lm_mp.c -lmpfr -lgmp
 * Usage:  euclid_lm_mp [--prec BITS] "a,b,c;d,e,f;..."  > out.obj
 *         euclid_lm_mp [--prec BITS] --batch --outdir DIR < netcodes.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpfr.h>

#ifndef MAXV
#define MAXV 2048   /* override: -DMAXV=... */
#endif
#define MAXF (2 * MAXV + 4)
#define MAXE (3 * MAXV + 6)
#define MAXDEG 16

static int NV, NF, NE;
static int FACES[MAXF][3];
static int EDGE_A[MAXE], EDGE_B[MAXE];
static int EDGE_IDX[MAXV + 1][MAXV + 1];
static int DIR_FACE[MAXV + 1][MAXV + 1];
static int VERT_DEG[MAXV + 1];
static int FLOWER_LEN[MAXV + 1];
static int FLOWER_E[MAXV + 1][MAXDEG];

/* specified flatfaces: frozen edges have bend identically 0 and are not
   LM variables. FREE_E lists the variable edges. */
static int IS_FLAT[MAXE];
static int FREE_E[MAXE];
static int NFREE;

static mpfr_prec_t PREC = 128;
static const char *FLATS_STR = NULL;   /* "a,b;c,d;..." vertex pairs */
static const char *SEED_PATH = NULL;   /* lines: a b bend-decimal */
static int NOGATE = 0;                 /* --nogate: LM without the dent gate */
static double ALPHA_DEG = 60.0;        /* --alpha: corner angle in degrees; 60 = Euclidean,
                                          smaller = hyperbolic, -> 0 = ideal */
static const char *DENTS_STR = NULL;   /* --dents "v1,v2": these must stay dented */
static int DENT_REQ[MAXV + 1];
static int NDENTREQ = 0;
static int DENTS_EXACT = 0;            /* --dents-exact: T<0 on listed, T>=0 off */
static int BENDS_ONLY = 0;             /* omit v/f lines from output */
static int PROVE = 0;                  /* solve, classify, refreeze, certify */
static int USE_SPARSE = 1;             /* fixed-pattern sparse Cholesky is the
                                          default; --dense restores the dense path
                                          (same pivot order, same operand order:
                                          the two are byte-identity-swept) */
static double PROVE_FLAT_TOL = 1e-8;
static int IS_PI[MAXE];                /* pancake case: edges pinned at gem/2 */
static int PANCAKE = 0;
static const char *NET_NAME = NULL;
static int CERT_OK;
static double CERT_SIG, CERT_H, CERT_RADIUS, CERT_DROP, CERT_ETA,
    CERT_TMIN, CERT_TMAX;
static const char *CERT_WHY = "";
static double CERT_DROP_AT;
static double PROVE_TOL_USED;
static mpfr_t BSAVE[MAXE];

/* ---------------- topology: identical to euclid_clean.c ---------------- */

static int parse_netcode(const char *s)
{
    NV = 0; NF = 0;
    const char *p = s;
    while (*p) {
        int a, b, c;
        if (sscanf(p, "%d,%d,%d", &a, &b, &c) != 3) return -1;
        if (NF >= MAXF || a < 1 || a > MAXV || b < 1 || b > MAXV || c < 1 || c > MAXV)
            return -1;
        FACES[NF][0] = a; FACES[NF][1] = b; FACES[NF][2] = c;
        NF++;
        if (a > NV) NV = a;
        if (b > NV) NV = b;
        if (c > NV) NV = c;
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }
    return NF;
}

static int build_topology(void)
{
    NE = 0;
    for (int i = 0; i <= MAXV; i++) {
        VERT_DEG[i] = 0;
        for (int j = 0; j <= MAXV; j++) { EDGE_IDX[i][j] = -1; DIR_FACE[i][j] = -1; }
    }
    for (int fi = 0; fi < NF; fi++) {
        int v[3] = { FACES[fi][0], FACES[fi][1], FACES[fi][2] };
        for (int k = 0; k < 3; k++) {
            int u = v[k], w = v[(k + 1) % 3];
            DIR_FACE[u][w] = fi;
            int lo = u < w ? u : w, hi = u < w ? w : u;
            if (EDGE_IDX[lo][hi] < 0) {
                if (NE >= MAXE) return -1;
                EDGE_IDX[lo][hi] = NE;
                EDGE_A[NE] = lo; EDGE_B[NE] = hi;
                VERT_DEG[lo]++; VERT_DEG[hi]++;
                NE++;
            }
        }
    }
    for (int v = 1; v <= NV; v++) {
        int start = -1;
        for (int fi = 0; fi < NF; fi++)
            if (FACES[fi][0] == v || FACES[fi][1] == v || FACES[fi][2] == v) { start = fi; break; }
        if (start < 0) return -1;
        int cur = start, k = 0;
        do {
            int a = FACES[cur][0], b = FACES[cur][1], c = FACES[cur][2];
            int third = (v == a) ? c : (v == b) ? a : b;
            int lo = v < third ? v : third, hi = v < third ? third : v;
            if (k >= MAXDEG) return -1;
            FLOWER_E[v][k++] = EDGE_IDX[lo][hi];
            cur = DIR_FACE[v][third];
        } while (cur != start);
        FLOWER_LEN[v] = k;
    }
    return 0;
}

/* -------- wish init: same LS system, plain double Gaussian solve -------- */

static int wish_init_double(double *wish)
{
    static double M[MAXV][MAXV + 1];
    static double g[MAXE];

    for (int e = 0; e < NE; e++)
        g[e] = 1.0 / VERT_DEG[EDGE_A[e]] + 1.0 / VERT_DEG[EDGE_B[e]];

    for (int i = 0; i < NV; i++)
        for (int j = 0; j <= NV; j++) M[i][j] = 0.0;
    for (int i = 1; i <= NV; i++) {
        M[i - 1][i - 1] = (double)VERT_DEG[i];
        for (int j = 1; j <= NV; j++) {
            if (i == j) continue;
            int lo = i < j ? i : j, hi = i < j ? j : i;
            if (EDGE_IDX[lo][hi] >= 0) M[i - 1][j - 1] = 1.0;
        }
    }
    /* per-vertex turning target: +gem, or -gem at --dents vertices, so
       the wish is the L2-min bend vector for the requested sign pattern
       (the undented pipeline is the all-plus special case) */
    for (int i = 0; i < NV; i++) M[i][NV] = DENT_REQ[i + 1] ? 2.0 : -2.0;
    for (int e = 0; e < NE; e++) {
        M[EDGE_A[e] - 1][NV] += g[e];
        M[EDGE_B[e] - 1][NV] += g[e];
    }
    for (int c = 0; c < NV; c++) {              /* partial-pivot elimination */
        int p = c;
        for (int r = c + 1; r < NV; r++)
            if (fabs(M[r][c]) > fabs(M[p][c])) p = r;
        if (fabs(M[p][c]) < 1e-14) return -1;
        if (p != c)
            for (int k = c; k <= NV; k++) { double t = M[c][k]; M[c][k] = M[p][k]; M[p][k] = t; }
        for (int r = 0; r < NV; r++) {
            if (r == c || M[r][c] == 0.0) continue;
            double f = M[r][c] / M[c][c];
            for (int k = c; k <= NV; k++) M[r][k] -= f * M[c][k];
        }
    }
    for (int e = 0; e < NE; e++) {
        double la = M[EDGE_A[e] - 1][NV] / M[EDGE_A[e] - 1][EDGE_A[e] - 1];
        double lb = M[EDGE_B[e] - 1][NV] / M[EDGE_B[e] - 1][EDGE_B[e] - 1];
        wish[e] = 2.0 * M_PI * 0.5 * (g[e] - la - lb);
    }
    return 0;
}

/* ---------------- MPFR working state ---------------- */

typedef mpfr_t Quat[4];

static mpfr_t BEND[MAXE], BEND_T[MAXE];
static mpfr_t R_BUF[3 * MAXV], RT_BUF[3 * MAXV];
static mpfr_t JCOL[MAXE][6];      /* column e: 3 rows at EDGE_A, 3 at EDGE_B */
static mpfr_t (*A_BUF)[MAXE];  /* J^T J (+ lambda D); heap: MAXE x MAXE
                                  static BSS would break the linker at
                                  MAXV 2048; pages are touched lazily */
static mpfr_t (*M_BUF)[MAXE];
static mpfr_t G_BUF[MAXE], D_BUF[MAXE], DELTA[MAXE];
static mpfr_t VXYZ[MAXV + 1][3];
static mpfr_t ALPHA, T1, T2, T3, T4, NORM, NTRIAL, LAMBDA, TOL;
static Quat QTMP1, QTMP2, QSTEP, QD;
static Quat QS[MAXDEG], DQS[MAXDEG], PPRE[MAXDEG + 1], SSUF[MAXDEG + 1];
static mpfr_t W1[3], W2[3], W3[3], W4[3], W5[3], W6[3];

/* big per-edge/per-vertex arrays are initialized lazily up to the
   sizes actually parsed (A_BUF/M_BUF are MAXE x MAXE static -- at
   MAXV 2048 initializing them in full would allocate GBs and take
   minutes; untouched static pages cost nothing). Monotone: batch
   inputs of growing size extend the initialized region. */
static int INITED_E = 0, INITED_V = 0;

static void mp_init_upto(int ne, int nv)
{
    if (ne > MAXE) ne = MAXE;
    if (nv > MAXV) nv = MAXV;
    if (!A_BUF && !USE_SPARSE) {
        A_BUF = calloc((size_t)MAXE * MAXE, sizeof(mpfr_t));
        M_BUF = calloc((size_t)MAXE * MAXE, sizeof(mpfr_t));
        if (!A_BUF || !M_BUF) { fprintf(stderr, "alloc failed\n"); exit(1); }
    }
    for (int e = INITED_E; e < ne; e++) {
        mpfr_init2(BEND[e], PREC); mpfr_init2(BEND_T[e], PREC);
        mpfr_init2(G_BUF[e], PREC); mpfr_init2(D_BUF[e], PREC); mpfr_init2(DELTA[e], PREC);
        for (int k = 0; k < 6; k++) mpfr_init2(JCOL[e][k], PREC);
    }
    /* A_BUF/M_BUF: new rows in full, then new columns of old rows */
    if (!USE_SPARSE) {
        for (int p = INITED_E; p < ne; p++)
            for (int q = 0; q < ne; q++) { mpfr_init2(A_BUF[p][q], PREC); mpfr_init2(M_BUF[p][q], PREC); }
        for (int p = 0; p < INITED_E; p++)
            for (int q = INITED_E; q < ne; q++) { mpfr_init2(A_BUF[p][q], PREC); mpfr_init2(M_BUF[p][q], PREC); }
    }
    for (int i = 3 * INITED_V; i < 3 * nv; i++) { mpfr_init2(R_BUF[i], PREC); mpfr_init2(RT_BUF[i], PREC); }
    for (int v = (INITED_V ? INITED_V + 1 : 0); v <= nv; v++)
        for (int k = 0; k < 3; k++) mpfr_init2(VXYZ[v][k], PREC);
    if (ne > INITED_E) INITED_E = ne;
    if (nv > INITED_V) INITED_V = nv;
}

static void mp_init_all(void)
{
    mpfr_inits2(PREC, ALPHA, T1, T2, T3, T4, NORM, NTRIAL, LAMBDA, TOL, (mpfr_ptr)0);
    for (int k = 0; k < 4; k++) {
        mpfr_init2(QTMP1[k], PREC); mpfr_init2(QTMP2[k], PREC);
        mpfr_init2(QSTEP[k], PREC); mpfr_init2(QD[k], PREC);
    }
    for (int t = 0; t < MAXDEG; t++)
        for (int k = 0; k < 4; k++) { mpfr_init2(QS[t][k], PREC); mpfr_init2(DQS[t][k], PREC); }
    for (int t = 0; t <= MAXDEG; t++)
        for (int k = 0; k < 4; k++) { mpfr_init2(PPRE[t][k], PREC); mpfr_init2(SSUF[t][k], PREC); }
    for (int k = 0; k < 3; k++)
        mpfr_inits2(PREC, W1[k], W2[k], W3[k], W4[k], W5[k], W6[k], (mpfr_ptr)0);
}

static void q_mul(const Quat a, const Quat b, Quat r)  /* r may alias a or b */
{
    static mpfr_t w, x, y, z, u;
    static int init = 0;
    if (!init) { mpfr_inits2(PREC, w, x, y, z, u, (mpfr_ptr)0); init = 1; }
    /* w = a0 b0 - a1 b1 - a2 b2 - a3 b3, etc. */
    mpfr_mul(w, a[0], b[0], MPFR_RNDN);
    mpfr_mul(u, a[1], b[1], MPFR_RNDN); mpfr_sub(w, w, u, MPFR_RNDN);
    mpfr_mul(u, a[2], b[2], MPFR_RNDN); mpfr_sub(w, w, u, MPFR_RNDN);
    mpfr_mul(u, a[3], b[3], MPFR_RNDN); mpfr_sub(w, w, u, MPFR_RNDN);
    mpfr_mul(x, a[0], b[1], MPFR_RNDN);
    mpfr_mul(u, a[1], b[0], MPFR_RNDN); mpfr_add(x, x, u, MPFR_RNDN);
    mpfr_mul(u, a[2], b[3], MPFR_RNDN); mpfr_add(x, x, u, MPFR_RNDN);
    mpfr_mul(u, a[3], b[2], MPFR_RNDN); mpfr_sub(x, x, u, MPFR_RNDN);
    mpfr_mul(y, a[0], b[2], MPFR_RNDN);
    mpfr_mul(u, a[1], b[3], MPFR_RNDN); mpfr_sub(y, y, u, MPFR_RNDN);
    mpfr_mul(u, a[2], b[0], MPFR_RNDN); mpfr_add(y, y, u, MPFR_RNDN);
    mpfr_mul(u, a[3], b[1], MPFR_RNDN); mpfr_add(y, y, u, MPFR_RNDN);
    mpfr_mul(z, a[0], b[3], MPFR_RNDN);
    mpfr_mul(u, a[1], b[2], MPFR_RNDN); mpfr_add(z, z, u, MPFR_RNDN);
    mpfr_mul(u, a[2], b[1], MPFR_RNDN); mpfr_sub(z, z, u, MPFR_RNDN);
    mpfr_mul(u, a[3], b[0], MPFR_RNDN); mpfr_add(z, z, u, MPFR_RNDN);
    mpfr_set(r[0], w, MPFR_RNDN); mpfr_set(r[1], x, MPFR_RNDN);
    mpfr_set(r[2], y, MPFR_RNDN); mpfr_set(r[3], z, MPFR_RNDN);
}

/* q_step(alpha, beta): (ca cb, -ca sb, -sa sb, sa cb), half angles */
static void q_step(const mpfr_t beta, Quat q, Quat dq, int want_d)
{
    static mpfr_t ca, sa, cb, sb, hb;
    static int init = 0;
    if (!init) { mpfr_inits2(PREC, ca, sa, cb, sb, hb, (mpfr_ptr)0); init = 1; }
    mpfr_div_ui(hb, ALPHA, 2, MPFR_RNDN);
    mpfr_cos(ca, hb, MPFR_RNDN); mpfr_sin(sa, hb, MPFR_RNDN);
    mpfr_div_ui(hb, beta, 2, MPFR_RNDN);
    mpfr_cos(cb, hb, MPFR_RNDN); mpfr_sin(sb, hb, MPFR_RNDN);
    mpfr_mul(q[0], ca, cb, MPFR_RNDN);
    mpfr_mul(q[1], ca, sb, MPFR_RNDN); mpfr_neg(q[1], q[1], MPFR_RNDN);
    mpfr_mul(q[2], sa, sb, MPFR_RNDN); mpfr_neg(q[2], q[2], MPFR_RNDN);
    mpfr_mul(q[3], sa, cb, MPFR_RNDN);
    if (want_d) {
        mpfr_mul(dq[0], ca, sb, MPFR_RNDN); mpfr_div_si(dq[0], dq[0], -2, MPFR_RNDN);
        mpfr_mul(dq[1], ca, cb, MPFR_RNDN); mpfr_div_si(dq[1], dq[1], -2, MPFR_RNDN);
        mpfr_mul(dq[2], sa, cb, MPFR_RNDN); mpfr_div_si(dq[2], dq[2], -2, MPFR_RNDN);
        mpfr_mul(dq[3], sa, sb, MPFR_RNDN); mpfr_div_si(dq[3], dq[3], -2, MPFR_RNDN);
    }
}

static void residual(mpfr_t *bend, mpfr_t *r)
{
    for (int v = 1; v <= NV; v++) {
        Quat Q;
        for (int k = 0; k < 4; k++) mpfr_set(QTMP1[k], k == 0 ? T1 : T2, MPFR_RNDN);
        mpfr_set_ui(QTMP1[0], 1, MPFR_RNDN);
        for (int k = 1; k < 4; k++) mpfr_set_ui(QTMP1[k], 0, MPFR_RNDN);
        (void)Q;
        int k = FLOWER_LEN[v];
        for (int t = 0; t < k; t++) {
            q_step(bend[FLOWER_E[v][t]], QSTEP, QD, 0);
            q_mul(QTMP1, QSTEP, QTMP1);
        }
        mpfr_set(r[3 * (v - 1) + 0], QTMP1[1], MPFR_RNDN);
        mpfr_set(r[3 * (v - 1) + 1], QTMP1[2], MPFR_RNDN);
        mpfr_set(r[3 * (v - 1) + 2], QTMP1[3], MPFR_RNDN);
    }
}

static void vec_norm2(mpfr_t out, mpfr_t *r, int n)
{
    mpfr_set_ui(out, 0, MPFR_RNDN);
    for (int i = 0; i < n; i++) {
        mpfr_mul(T1, r[i], r[i], MPFR_RNDN);
        mpfr_add(out, out, T1, MPFR_RNDN);
    }
    mpfr_sqrt(out, out, MPFR_RNDN);
}

/* Jacobian: per vertex, dQ = P[t] * dq[t] * S[t+1]; column e gets 3 rows
   at each of its two endpoint vertices. Stored in JCOL[e][0..5]:
   rows 3*(EDGE_A[e]-1)+i and 3*(EDGE_B[e]-1)+i. */
static void jacobian_cols(mpfr_t *bend)
{
    for (int e = 0; e < NE; e++)
        for (int k = 0; k < 6; k++) mpfr_set_ui(JCOL[e][k], 0, MPFR_RNDN);

    for (int v = 1; v <= NV; v++) {
        int k = FLOWER_LEN[v];
        for (int t = 0; t < k; t++)
            q_step(bend[FLOWER_E[v][t]], QS[t], DQS[t], 1);
        mpfr_set_ui(PPRE[0][0], 1, MPFR_RNDN);
        for (int c = 1; c < 4; c++) mpfr_set_ui(PPRE[0][c], 0, MPFR_RNDN);
        for (int t = 0; t < k; t++) q_mul(PPRE[t], QS[t], PPRE[t + 1]);
        mpfr_set_ui(SSUF[k][0], 1, MPFR_RNDN);
        for (int c = 1; c < 4; c++) mpfr_set_ui(SSUF[k][c], 0, MPFR_RNDN);
        for (int t = k - 1; t >= 0; t--) q_mul(QS[t], SSUF[t + 1], SSUF[t]);

        for (int t = 0; t < k; t++) {
            int e = FLOWER_E[v][t];
            if (IS_FLAT[e]) continue;          /* frozen: not a variable */
            q_mul(PPRE[t], DQS[t], QTMP2);
            q_mul(QTMP2, SSUF[t + 1], QTMP2);
            int side = (EDGE_A[e] == v) ? 0 : 3;
            for (int c = 0; c < 3; c++)
                mpfr_add(JCOL[e][side + c], JCOL[e][side + c], QTMP2[c + 1], MPFR_RNDN);
        }
    }
}

/* ================= sparse fixed-pattern Cholesky path =================
 * --sparse: same pivot order (natural), same per-entry k-accumulation
 * order as the dense chol_solve_neg below.  Every operation skipped is
 * one whose dense operands are structural zeros, which in mpfr are
 * exact 0.0 (mul by exact 0 gives exact 0; subtracting exact 0 is
 * exact), so no computed value can differ: the factor, the solves, and
 * hence the LM iterates are bit-identical to the dense path.  The only
 * change that would break bit-identity is pivot reordering, which this
 * stage does not do.  Pattern = free-edge adjacency (columns overlap
 * iff edges share an endpoint; see normal_matrix) plus symbolic fill
 * in natural order; rebuilt whenever NFREE/FREE_E change (freeze
 * ladder rungs). */
static int SP_N = -1;                     /* pattern built for this NFREE */
static int SP_FREE_SIG[MAXE];             /* FREE_E copy at build time */
static int *SPR_PTR, *SPR_COL;            /* strict lower rows: i -> cols k<i asc */
static int *SPC_PTR, *SPC_ROW, *SPC_SLOT; /* strict lower cols: j -> rows i>j asc */
static mpfr_t *SPA_V, *SPM_V;             /* values on the row slots */
static mpfr_t *SPA_D, *SPM_D;             /* diagonals */
static int SP_NNZ = 0, SP_MP_SLOTS = 0, SP_MP_DIAG = 0;

static void sparse_build_pattern(void)
{
    int n = NFREE;
    int words = (n + 63) / 64;
    unsigned long long *col = calloc((size_t)n * words, sizeof(*col));
    int *efree = malloc(sizeof(int) * (MAXE + 1));
    if (!col || !efree) { fprintf(stderr, "sparse alloc failed\n"); exit(1); }
    for (int e = 0; e < NE; e++) efree[e] = -1;
    for (int pf = 0; pf < n; pf++) efree[FREE_E[pf]] = pf;
    /* A-pattern: free edges sharing a vertex */
    for (int v = 1; v <= NV; v++)
        for (int s = 0; s < FLOWER_LEN[v]; s++) {
            int pf = efree[FLOWER_E[v][s]];
            if (pf < 0) continue;
            for (int t = s + 1; t < FLOWER_LEN[v]; t++) {
                int qf = efree[FLOWER_E[v][t]];
                if (qf < 0) continue;
                int lo = pf < qf ? pf : qf, hi = pf < qf ? qf : pf;
                col[(size_t)lo * words + hi / 64] |= 1ULL << (hi % 64);
            }
        }
    /* symbolic fill, natural order */
    for (int k = 0; k < n; k++) {
        unsigned long long *ck = col + (size_t)k * words;
        int p = -1;
        for (int w = k / 64; w < words && p < 0; w++)
            if (ck[w]) {
                unsigned long long m = ck[w];
                if (w == k / 64) m &= ~((k % 64 == 63) ? 0ULL : ((1ULL << (k % 64 + 1)) - 1));
                if (m) p = w * 64 + __builtin_ctzll(m);
            }
        if (p < 0) continue;
        unsigned long long *cp = col + (size_t)p * words;
        for (int w = 0; w < words; w++) cp[w] |= ck[w];
        cp[p / 64] &= ~(1ULL << (p % 64));
    }
    /* row/col lists (append in k-ascending scan => sorted) */
    free(SPR_PTR); free(SPR_COL); free(SPC_PTR); free(SPC_ROW); free(SPC_SLOT);
    int *rcnt = calloc(n + 1, sizeof(int));
    int nnz = 0;
    for (int k = 0; k < n; k++) {
        unsigned long long *ck = col + (size_t)k * words;
        for (int w = 0; w < words; w++)
            for (unsigned long long m = ck[w]; m; m &= m - 1) {
                int i = w * 64 + __builtin_ctzll(m);
                if (i > k) { rcnt[i]++; nnz++; }
            }
    }
    SPR_PTR = malloc(sizeof(int) * (n + 1)); SPR_COL = malloc(sizeof(int) * (nnz ? nnz : 1));
    SPC_PTR = malloc(sizeof(int) * (n + 1)); SPC_ROW = malloc(sizeof(int) * (nnz ? nnz : 1));
    SPC_SLOT = malloc(sizeof(int) * (nnz ? nnz : 1));
    if (!SPR_PTR || !SPR_COL || !SPC_PTR || !SPC_ROW || !SPC_SLOT) { fprintf(stderr, "sparse alloc failed\n"); exit(1); }
    SPR_PTR[0] = 0;
    for (int i = 0; i < n; i++) SPR_PTR[i + 1] = SPR_PTR[i] + rcnt[i];
    int *rfill = calloc(n, sizeof(int));
    int *ccnt = calloc(n + 1, sizeof(int));
    for (int k = 0; k < n; k++) {                /* k asc => each row's cols asc */
        unsigned long long *ck = col + (size_t)k * words;
        for (int w = 0; w < words; w++)
            for (unsigned long long m = ck[w]; m; m &= m - 1) {
                int i = w * 64 + __builtin_ctzll(m);
                if (i > k) { SPR_COL[SPR_PTR[i] + rfill[i]++] = k; ccnt[k]++; }
            }
    }
    SPC_PTR[0] = 0;
    for (int k = 0; k < n; k++) SPC_PTR[k + 1] = SPC_PTR[k] + ccnt[k];
    int *cfill = calloc(n, sizeof(int));
    for (int i = 0; i < n; i++)                  /* i asc => each col's rows asc */
        for (int s = SPR_PTR[i]; s < SPR_PTR[i + 1]; s++) {
            int k = SPR_COL[s];
            int c = SPC_PTR[k] + cfill[k]++;
            SPC_ROW[c] = i; SPC_SLOT[c] = s;
        }
    free(col); free(efree); free(rcnt); free(rfill); free(ccnt); free(cfill);
    /* mpfr slots: grow-only */
    if (!SPA_V) {
        SPA_V = calloc(MAXE * 32, sizeof(mpfr_t)); SPM_V = calloc(MAXE * 32, sizeof(mpfr_t));
        SPA_D = calloc(MAXE, sizeof(mpfr_t)); SPM_D = calloc(MAXE, sizeof(mpfr_t));
        if (!SPA_V || !SPM_V || !SPA_D || !SPM_D) { fprintf(stderr, "sparse alloc failed\n"); exit(1); }
    }
    if (nnz > MAXE * 32) { fprintf(stderr, "sparse fill %d exceeds slot cap %d\n", nnz, MAXE * 32); exit(1); }
    for (int s = SP_MP_SLOTS; s < nnz; s++) { mpfr_init2(SPA_V[s], PREC); mpfr_init2(SPM_V[s], PREC); }
    if (nnz > SP_MP_SLOTS) SP_MP_SLOTS = nnz;
    for (int i = SP_MP_DIAG; i < n; i++) { mpfr_init2(SPA_D[i], PREC); mpfr_init2(SPM_D[i], PREC); }
    if (n > SP_MP_DIAG) SP_MP_DIAG = n;
    SP_NNZ = nnz;
    SP_N = n;
    memcpy(SP_FREE_SIG, FREE_E, sizeof(int) * n);
}

static void sparse_check_pattern(void)
{
    if (SP_N == NFREE && memcmp(SP_FREE_SIG, FREE_E, sizeof(int) * NFREE) == 0)
        return;
    sparse_build_pattern();
}

/* same per-entry operand order as normal_matrix: entry (i,j), i>j, is
   the dense upper-triangle entry A[pf=j][qf=i] (then mirrored), so the
   sp/sq endpoint loops run with p = FREE_E[j], q = FREE_E[i]. */
static void sparse_entry(mpfr_t dst, int p, int q)
{
    mpfr_set_ui(dst, 0, MPFR_RNDN);
    for (int sp = 0; sp < 2; sp++) {
        int vp = sp ? EDGE_B[p] : EDGE_A[p];
        for (int sq = 0; sq < 2; sq++) {
            int vq = sq ? EDGE_B[q] : EDGE_A[q];
            if (vp != vq) continue;
            for (int c = 0; c < 3; c++) {
                mpfr_mul(T1, JCOL[p][3 * sp + c], JCOL[q][3 * sq + c], MPFR_RNDN);
                mpfr_add(dst, dst, T1, MPFR_RNDN);
            }
        }
    }
}

static void normal_matrix_sparse(void)
{
    sparse_check_pattern();
    for (int i = 0; i < SP_N; i++) {
        sparse_entry(SPA_D[i], FREE_E[i], FREE_E[i]);
        for (int s = SPR_PTR[i]; s < SPR_PTR[i + 1]; s++)
            sparse_entry(SPA_V[s], FREE_E[SPR_COL[s]], FREE_E[i]);
    }
}

static int chol_solve_neg_sparse(void)
{
    int n = SP_N;
    for (int i = 0; i < n; i++) {
        mpfr_set(SPM_D[i], SPA_D[i], MPFR_RNDN);
        mpfr_mul(T1, LAMBDA, D_BUF[i], MPFR_RNDN);
        mpfr_add(SPM_D[i], SPM_D[i], T1, MPFR_RNDN);
    }
    for (int s = 0; s < SP_NNZ; s++) mpfr_set(SPM_V[s], SPA_V[s], MPFR_RNDN);
    for (int j = 0; j < n; j++) {
        for (int s = SPR_PTR[j]; s < SPR_PTR[j + 1]; s++) {
            mpfr_mul(T1, SPM_V[s], SPM_V[s], MPFR_RNDN);
            mpfr_sub(SPM_D[j], SPM_D[j], T1, MPFR_RNDN);
        }
        if (mpfr_sgn(SPM_D[j]) <= 0) return -1;
        mpfr_sqrt(SPM_D[j], SPM_D[j], MPFR_RNDN);
        for (int c = SPC_PTR[j]; c < SPC_PTR[j + 1]; c++) {
            int i = SPC_ROW[c], sij = SPC_SLOT[c];
            int a = SPR_PTR[i], ae = SPR_PTR[i + 1];
            int b = SPR_PTR[j], be = SPR_PTR[j + 1];
            while (a < ae && b < be) {           /* merge: k ascending */
                int ka = SPR_COL[a], kb = SPR_COL[b];
                if (ka == kb) {
                    mpfr_mul(T1, SPM_V[a], SPM_V[b], MPFR_RNDN);
                    mpfr_sub(SPM_V[sij], SPM_V[sij], T1, MPFR_RNDN);
                    a++; b++;
                } else if (ka < kb) a++;
                else b++;
            }
            mpfr_div(SPM_V[sij], SPM_V[sij], SPM_D[j], MPFR_RNDN);
        }
    }
    for (int i = 0; i < n; i++) {                /* forward: L y = -g */
        mpfr_neg(DELTA[i], G_BUF[i], MPFR_RNDN);
        for (int s = SPR_PTR[i]; s < SPR_PTR[i + 1]; s++) {
            mpfr_mul(T1, SPM_V[s], DELTA[SPR_COL[s]], MPFR_RNDN);
            mpfr_sub(DELTA[i], DELTA[i], T1, MPFR_RNDN);
        }
        mpfr_div(DELTA[i], DELTA[i], SPM_D[i], MPFR_RNDN);
    }
    for (int i = n - 1; i >= 0; i--) {           /* back: L^T x = y */
        for (int c = SPC_PTR[i]; c < SPC_PTR[i + 1]; c++) {
            mpfr_mul(T1, SPM_V[SPC_SLOT[c]], DELTA[SPC_ROW[c]], MPFR_RNDN);
            mpfr_sub(DELTA[i], DELTA[i], T1, MPFR_RNDN);
        }
        mpfr_div(DELTA[i], DELTA[i], SPM_D[i], MPFR_RNDN);
    }
    return 0;
}
/* =============== end sparse fixed-pattern Cholesky path =============== */

/* normal matrix A = J^T J via shared rows: columns p, q overlap only on
   rows of shared endpoint vertices. */
static void normal_matrix(void)
{
    if (USE_SPARSE) { normal_matrix_sparse(); return; }
    for (int pf = 0; pf < NFREE; pf++)
        for (int qf = pf; qf < NFREE; qf++) mpfr_set_ui(A_BUF[pf][qf], 0, MPFR_RNDN);

    for (int pf = 0; pf < NFREE; pf++) {
        int p = FREE_E[pf];
        for (int qf = pf; qf < NFREE; qf++) {
            int q = FREE_E[qf];
            for (int sp = 0; sp < 2; sp++) {
                int vp = sp ? EDGE_B[p] : EDGE_A[p];
                for (int sq = 0; sq < 2; sq++) {
                    int vq = sq ? EDGE_B[q] : EDGE_A[q];
                    if (vp != vq) continue;
                    for (int c = 0; c < 3; c++) {
                        mpfr_mul(T1, JCOL[p][3 * sp + c], JCOL[q][3 * sq + c], MPFR_RNDN);
                        mpfr_add(A_BUF[pf][qf], A_BUF[pf][qf], T1, MPFR_RNDN);
                    }
                }
            }
        }
    }
    for (int pf = 0; pf < NFREE; pf++)
        for (int qf = 0; qf < pf; qf++) mpfr_set(A_BUF[pf][qf], A_BUF[qf][pf], MPFR_RNDN);
}

static void grad_vec(mpfr_t *r)
{
    for (int pf = 0; pf < NFREE; pf++) {
        int p = FREE_E[pf];
        mpfr_set_ui(G_BUF[pf], 0, MPFR_RNDN);
        for (int sp = 0; sp < 2; sp++) {
            int vp = sp ? EDGE_B[p] : EDGE_A[p];
            for (int c = 0; c < 3; c++) {
                mpfr_mul(T1, JCOL[p][3 * sp + c], r[3 * (vp - 1) + c], MPFR_RNDN);
                mpfr_add(G_BUF[pf], G_BUF[pf], T1, MPFR_RNDN);
            }
        }
    }
}

/* Cholesky solve of M x = -g, M = A + lambda*diag(D). Returns -1 if a
   pivot is nonpositive (raise lambda, like dgesv info != 0). */
static int chol_solve_neg(void)
{
    if (USE_SPARSE) return chol_solve_neg_sparse();
    int n = NFREE;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) mpfr_set(M_BUF[j][i], A_BUF[j][i], MPFR_RNDN);
        mpfr_mul(T1, LAMBDA, D_BUF[i], MPFR_RNDN);
        mpfr_add(M_BUF[i][i], M_BUF[i][i], T1, MPFR_RNDN);
    }
    for (int j = 0; j < n; j++) {                 /* lower Cholesky in place */
        for (int k = 0; k < j; k++) {
            mpfr_mul(T1, M_BUF[j][k], M_BUF[j][k], MPFR_RNDN);
            mpfr_sub(M_BUF[j][j], M_BUF[j][j], T1, MPFR_RNDN);
        }
        if (mpfr_sgn(M_BUF[j][j]) <= 0) return -1;
        mpfr_sqrt(M_BUF[j][j], M_BUF[j][j], MPFR_RNDN);
        for (int i = j + 1; i < n; i++) {
            for (int k = 0; k < j; k++) {
                mpfr_mul(T1, M_BUF[i][k], M_BUF[j][k], MPFR_RNDN);
                mpfr_sub(M_BUF[i][j], M_BUF[i][j], T1, MPFR_RNDN);
            }
            mpfr_div(M_BUF[i][j], M_BUF[i][j], M_BUF[j][j], MPFR_RNDN);
        }
    }
    for (int i = 0; i < n; i++) {                 /* forward: L y = -g */
        mpfr_neg(DELTA[i], G_BUF[i], MPFR_RNDN);
        for (int k = 0; k < i; k++) {
            mpfr_mul(T1, M_BUF[i][k], DELTA[k], MPFR_RNDN);
            mpfr_sub(DELTA[i], DELTA[i], T1, MPFR_RNDN);
        }
        mpfr_div(DELTA[i], DELTA[i], M_BUF[i][i], MPFR_RNDN);
    }
    for (int i = n - 1; i >= 0; i--) {            /* back: L^T x = y */
        for (int k = i + 1; k < n; k++) {
            mpfr_mul(T1, M_BUF[k][i], DELTA[k], MPFR_RNDN);
            mpfr_sub(DELTA[i], DELTA[i], T1, MPFR_RNDN);
        }
        mpfr_div(DELTA[i], DELTA[i], M_BUF[i][i], MPFR_RNDN);
    }
    return 0;
}

static int has_dent(mpfr_t *bend)
{
    for (int v = 1; v <= NV; v++) {
        mpfr_set_ui(T2, 0, MPFR_RNDN);
        for (int t = 0; t < FLOWER_LEN[v]; t++)
            mpfr_add(T2, T2, bend[FLOWER_E[v][t]], MPFR_RNDN);
        if (mpfr_sgn(T2) < 0) return v;
    }
    return 0;
}

/* --dents gate: every listed vertex must have negative turning; other
   vertices are unconstrained -- unless DENTS_EXACT, in which case
   unlisted vertices must have T >= 0 (the full sign pattern).
   Returns a violating vertex or 0. */
static int dents_gate_bad(mpfr_t *bend)
{
    if (DENTS_EXACT) {
        /* principal branch: no fold past a halfturn (|b| < pi); a bend
           beyond shut is never embedded and lets the unwrapped turning
           diverge from the record's */
        mpfr_const_pi(T2, MPFR_RNDN);
        for (int e = 0; e < NE; e++)
            if (mpfr_cmpabs(bend[e], T2) >= 0) return NV + 1 + e;
    }
    for (int v = 1; v <= NV; v++) {
        mpfr_set_ui(T2, 0, MPFR_RNDN);
        for (int t = 0; t < FLOWER_LEN[v]; t++)
            mpfr_add(T2, T2, bend[FLOWER_E[v][t]], MPFR_RNDN);
        /* closed orthant: listed vertices T <= 0, unlisted T >= 0; an
           exactly flat vertex passes both signs (dents of flat vertices
           are in the eye of the beholder) */
        if (DENT_REQ[v]) { if (mpfr_sgn(T2) > 0) return v; }
        else if (DENTS_EXACT) { if (mpfr_sgn(T2) < 0) return v; }
    }
    return 0;
}

static int parse_dents(const char *str)
{
    for (int v = 0; v <= MAXV; v++) DENT_REQ[v] = 0;
    NDENTREQ = 0;
    if (!str) return 0;
    const char *p = str;
    while (*p) {
        int v;
        if (sscanf(p, "%d", &v) != 1) return -1;
        if (v < 1 || v > NV) return -1;
        DENT_REQ[v] = 1; NDENTREQ++;
        while (*p && *p != ',') p++;
        if (*p == ',') p++;
    }
    return 0;
}

static int LM_ITERS = 0;

static int lm_solve(int maxiter)
{
    /* diagnostic veto counters (stderr only; stdout unchanged) */
    int dg_gate = 0, dg_gate_only = 0, dg_norm = 0, dg_chol = 0, dg_lastv = 0;
    int rows = 3 * NV;
    mpfr_set_d(LAMBDA, 1.0, MPFR_RNDN);

    residual(BEND, R_BUF);
    vec_norm2(NORM, R_BUF, rows);

    for (int it = 0; it < maxiter; it++) {
        LM_ITERS = it;
        if (mpfr_cmp(NORM, TOL) <= 0) return 0;

        jacobian_cols(BEND);
        normal_matrix();
        grad_vec(R_BUF);
        for (int pf = 0; pf < NFREE; pf++)
            mpfr_set(D_BUF[pf], USE_SPARSE ? SPA_D[pf] : A_BUF[pf][pf], MPFR_RNDN);

        int accepted = 0;
        for (int retry = 0; retry < 20; retry++) {
            if (chol_solve_neg() != 0) {
                dg_chol++;
                mpfr_mul_ui(LAMBDA, LAMBDA, 10, MPFR_RNDN);
                if (mpfr_cmp_d(LAMBDA, 1e12) > 0) goto diag_fail;
                continue;
            }
            for (int e = 0; e < NE; e++) mpfr_set(BEND_T[e], BEND[e], MPFR_RNDN);
            for (int pf = 0; pf < NFREE; pf++)
                mpfr_add(BEND_T[FREE_E[pf]], BEND[FREE_E[pf]], DELTA[pf], MPFR_RNDN);
            int gate_bad;
            if (NOGATE)                        gate_bad = 0;
            else if (NDENTREQ || DENTS_EXACT)  gate_bad = dents_gate_bad(BEND_T);
            else                               gate_bad = has_dent(BEND_T);
            residual(BEND_T, RT_BUF);
            vec_norm2(NTRIAL, RT_BUF, rows);

            if (gate_bad) { dg_gate++; dg_lastv = gate_bad;
                            if (mpfr_cmp(NTRIAL, NORM) < 0) dg_gate_only++; }
            else if (mpfr_cmp(NTRIAL, NORM) >= 0) dg_norm++;
            if (mpfr_cmp(NTRIAL, NORM) < 0 && gate_bad == 0) {
                for (int e = 0; e < NE; e++) mpfr_set(BEND[e], BEND_T[e], MPFR_RNDN);
                for (int i = 0; i < rows; i++) mpfr_set(R_BUF[i], RT_BUF[i], MPFR_RNDN);
                mpfr_set(NORM, NTRIAL, MPFR_RNDN);
                mpfr_mul_d(LAMBDA, LAMBDA, 0.3, MPFR_RNDN);
                if (mpfr_cmp_d(LAMBDA, 1e-30) < 0) mpfr_set_d(LAMBDA, 1e-30, MPFR_RNDN);
                accepted = 1;
                break;
            }
            mpfr_mul_ui(LAMBDA, LAMBDA, 10, MPFR_RNDN);
            if (mpfr_cmp_d(LAMBDA, 1e12) > 0) goto diag_fail;
        }
        if (!accepted) goto diag_fail;
    }
    if (mpfr_cmp(NORM, TOL) <= 0) return 0;
diag_fail:
    fprintf(stderr, "# lm_diag fail: it=%d norm=%.3e lambda=%.3e "
            "gate_vetoes=%d (would_have_accepted=%d) norm_vetoes=%d "
            "chol_fails=%d last_gate_vertex=%d\n",
            LM_ITERS, mpfr_get_d(NORM, MPFR_RNDN), mpfr_get_d(LAMBDA, MPFR_RNDN),
            dg_gate, dg_gate_only, dg_norm, dg_chol, dg_lastv);
    return -1;
}

/* ---------------- realize (MPFR) ---------------- */

static void v_sub(mpfr_t a[3], mpfr_t b[3], mpfr_t r[3])
{ for (int k = 0; k < 3; k++) mpfr_sub(r[k], a[k], b[k], MPFR_RNDN); }
static void v_addto(mpfr_t r[3], mpfr_t a[3])
{ for (int k = 0; k < 3; k++) mpfr_add(r[k], r[k], a[k], MPFR_RNDN); }
static void v_scale(const mpfr_t s, mpfr_t a[3], mpfr_t r[3])
{ for (int k = 0; k < 3; k++) mpfr_mul(r[k], a[k], s, MPFR_RNDN); }

static void place_third(int a, int b, int p, const mpfr_t theta, int c)
{
    /* m = (A+B)/2; u = (P-m)/|P-m|; v = u x e / |u x e|, e = B-A;
       C = m + (sqrt3/2)(-cos th * u + sin th * v) */
    for (int k = 0; k < 3; k++) {
        mpfr_add(W1[k], VXYZ[a][k], VXYZ[b][k], MPFR_RNDN);
        mpfr_div_ui(W1[k], W1[k], 2, MPFR_RNDN);          /* W1 = m */
    }
    v_sub(VXYZ[b], VXYZ[a], W2);                          /* W2 = e */
    v_sub(VXYZ[p], W1, W3);                               /* W3 = p_perp */
    mpfr_set_ui(T1, 0, MPFR_RNDN);
    for (int k = 0; k < 3; k++) {
        mpfr_mul(T2, W3[k], W3[k], MPFR_RNDN);
        mpfr_add(T1, T1, T2, MPFR_RNDN);
    }
    mpfr_sqrt(T1, T1, MPFR_RNDN);
    for (int k = 0; k < 3; k++) mpfr_div(W3[k], W3[k], T1, MPFR_RNDN);  /* u */
    /* W4 = u x e */
    mpfr_mul(W4[0], W3[1], W2[2], MPFR_RNDN); mpfr_mul(T1, W3[2], W2[1], MPFR_RNDN);
    mpfr_sub(W4[0], W4[0], T1, MPFR_RNDN);
    mpfr_mul(W4[1], W3[2], W2[0], MPFR_RNDN); mpfr_mul(T1, W3[0], W2[2], MPFR_RNDN);
    mpfr_sub(W4[1], W4[1], T1, MPFR_RNDN);
    mpfr_mul(W4[2], W3[0], W2[1], MPFR_RNDN); mpfr_mul(T1, W3[1], W2[0], MPFR_RNDN);
    mpfr_sub(W4[2], W4[2], T1, MPFR_RNDN);
    mpfr_set_ui(T1, 0, MPFR_RNDN);
    for (int k = 0; k < 3; k++) {
        mpfr_mul(T2, W4[k], W4[k], MPFR_RNDN);
        mpfr_add(T1, T1, T2, MPFR_RNDN);
    }
    mpfr_sqrt(T1, T1, MPFR_RNDN);
    for (int k = 0; k < 3; k++) mpfr_div(W4[k], W4[k], T1, MPFR_RNDN);  /* v */
    mpfr_cos(T2, theta, MPFR_RNDN); mpfr_neg(T2, T2, MPFR_RNDN);
    mpfr_sin(T3, theta, MPFR_RNDN);
    v_scale(T2, W3, W5);
    v_scale(T3, W4, W6);
    v_addto(W5, W6);
    mpfr_sqrt_ui(T4, 3, MPFR_RNDN); mpfr_div_ui(T4, T4, 2, MPFR_RNDN);
    v_scale(T4, W5, W5);
    for (int k = 0; k < 3; k++)
        mpfr_add(VXYZ[c][k], W1[k], W5[k], MPFR_RNDN);
}

static int realize(void)
{
    static int queue[MAXF], placed_f[MAXF];
    int b0 = FACES[0][0], b1 = FACES[0][1], b2 = FACES[0][2];
    for (int fi = 0; fi < NF; fi++) placed_f[fi] = 0;
    for (int k = 0; k < 3; k++) {
        mpfr_set_ui(VXYZ[b0][k], 0, MPFR_RNDN);
        mpfr_set_ui(VXYZ[b1][k], 0, MPFR_RNDN);
        mpfr_set_ui(VXYZ[b2][k], 0, MPFR_RNDN);
    }
    mpfr_set_d(VXYZ[b0][2], 0.5, MPFR_RNDN);
    mpfr_set_d(VXYZ[b1][2], -0.5, MPFR_RNDN);
    mpfr_sqrt_ui(T1, 3, MPFR_RNDN); mpfr_div_ui(T1, T1, 2, MPFR_RNDN);
    mpfr_set(VXYZ[b2][0], T1, MPFR_RNDN);
    placed_f[0] = 1;
    int qh = 0, qt = 0;
    queue[qt++] = 0;
    while (qh < qt) {
        int fi = queue[qh++];
        int v[3] = { FACES[fi][0], FACES[fi][1], FACES[fi][2] };
        for (int k = 0; k < 3; k++) {
            int a = v[k], b = v[(k + 1) % 3];
            int other_fi = DIR_FACE[b][a];
            if (placed_f[other_fi]) continue;
            int oa = FACES[other_fi][0], ob = FACES[other_fi][1], oc = FACES[other_fi][2];
            int c = (oa != a && oa != b) ? oa : (ob != a && ob != b) ? ob : oc;
            int p = (v[0] != a && v[0] != b) ? v[0] : (v[1] != a && v[1] != b) ? v[1] : v[2];
            int lo = a < b ? a : b, hi = a < b ? b : a;
            place_third(a, b, p, BEND[EDGE_IDX[lo][hi]], c);
            placed_f[other_fi] = 1;
            queue[qt++] = other_fi;
        }
    }
    return 0;
}

static int write_obj(FILE *fp)
{
    mpfr_fprintf(fp, "# euclid_lm_mp prec=%ld resid=%.3Re iters=%d nflat=%d\n",
                 (long)PREC, NORM, LM_ITERS, NE - NFREE);
    if (PROVE)
        fprintf(fp, "# cert ok %d sig_low %.9e h %.6e radius %.6e drop %.6e "
                    "eta %.6e tmin %.9e tmax %.9e why %s\n",
                CERT_OK, CERT_SIG, CERT_H, CERT_RADIUS, CERT_DROP,
                CERT_ETA, CERT_TMIN, CERT_TMAX,
                CERT_OK ? "-" : CERT_WHY);
    for (int e = 0; e < NE; e++) {
        if (IS_FLAT[e])
            fprintf(fp, "# bend %d %d %d 0 flat\n", e, EDGE_A[e], EDGE_B[e]);
        else
            mpfr_fprintf(fp, "# bend %d %d %d %.40Re\n", e, EDGE_A[e], EDGE_B[e], BEND[e]);
    }
    if (BENDS_ONLY) return ferror(fp) ? -1 : 0;
    for (int v = 1; v <= NV; v++)
        mpfr_fprintf(fp, "v %.40Rf %.40Rf %.40Rf\n",
                     VXYZ[v][0], VXYZ[v][1], VXYZ[v][2]);
    for (int fi = 0; fi < NF; fi++)
        fprintf(fp, "f %d %d %d\n", FACES[fi][0], FACES[fi][1], FACES[fi][2]);
    return ferror(fp) ? -1 : 0;
}

/* ============================================================
 *   Take-1 certificate (mechanical port of bendprover/mpkernel.py).
 *
 *   Pair variables per fold edge e: (c_e, s_e) initialized to
 *   (cos b_e/2, -sin b_e/2), so the vertex word
 *       prod over star of T60 (x) (c_e, s_e, 0, 0)   =  eps_v * 1
 *   is exactly the word the LM solver drove to zero (its q_step is
 *   q_z(alpha) q_x(-beta)). Unit rows c^2 + s^2 = 1. Newton polish,
 *   then: numerically selected square row subsystem, Cholesky factor
 *   witness sigma lower bound, Kantorovich with L = 4*arity^2,
 *   rounding slop SLOP * ulp per the take-1 model. Radius is in pair
 *   space; |bend error| <= 2*radius/(1-radius) for radius < 1.
 * ============================================================ */

#define MAXP (2 * MAXE)
#define NRMAX (4 * MAXV + MAXE)
#define CSLOP 64.0

static mpfr_t PC[MAXE], PS[MAXE];          /* pair values per fold edge */
static int SGN[MAXV + 1];                  /* eps_v */
static int ACT[MAXV], NACT;                /* active vertices */
static int FIDX[MAXE];                     /* edge -> fold index, -1 */
static mpfr_t R2[NRMAX], G2[MAXP], DX2[MAXP];
static mpfr_t JB[MAXE][2][4][2];           /* fold e, endpoint side, comp, d{c,s} */
static mpfr_t (*A2)[MAXP], (*C2)[MAXP], (*B2)[MAXP];  /* heap: MAXP^2
   statics would put ~15GB in BSS at MAXV 2048 (ADRP range error);
   allocated on demand, pages touched lazily */
static double (*JD)[MAXP];
static int CERT_INITED_P = 0;

static void cert_init_upto(int np)
{
    if (np > MAXP) np = MAXP;
    if (!A2) {
        A2 = calloc((size_t)MAXP * MAXP, sizeof(mpfr_t));
        C2 = calloc((size_t)MAXP * MAXP, sizeof(mpfr_t));
        B2 = calloc((size_t)MAXP * MAXP, sizeof(mpfr_t));
        JD = calloc((size_t)NRMAX * MAXP, sizeof(double));
        if (!A2 || !C2 || !B2 || !JD) { fprintf(stderr, "cert alloc failed\n"); exit(1); }
    }
    for (int i = CERT_INITED_P; i < np; i++)
        for (int j = 0; j < np; j++) {
            mpfr_init2(A2[i][j], PREC); mpfr_init2(C2[i][j], PREC); mpfr_init2(B2[i][j], PREC);
        }
    for (int i = 0; i < CERT_INITED_P; i++)
        for (int j = CERT_INITED_P; j < np; j++) {
            mpfr_init2(A2[i][j], PREC); mpfr_init2(C2[i][j], PREC); mpfr_init2(B2[i][j], PREC);
        }
    if (np > CERT_INITED_P) CERT_INITED_P = np;
}
static int SEL[MAXP], INSEL[NRMAX];
static mpfr_t T60Q[4], QA[4], QB[4], QC[4];
static mpfr_t PRE[2 * MAXDEG + 2][4], SUF[2 * MAXDEG + 2][4];
static mpfr_t CT1, CT2, CT3;

static void cert_init_all(void)
{
    for (int e = 0; e < MAXE; e++) {
        mpfr_init2(PC[e], PREC); mpfr_init2(PS[e], PREC);
        for (int s = 0; s < 2; s++)
            for (int c = 0; c < 4; c++)
                for (int k = 0; k < 2; k++) mpfr_init2(JB[e][s][c][k], PREC);
    }
    for (int i = 0; i < NRMAX; i++) mpfr_init2(R2[i], PREC);
    for (int i = 0; i < MAXP; i++) { mpfr_init2(G2[i], PREC); mpfr_init2(DX2[i], PREC); }
    for (int k = 0; k < 4; k++) {
        mpfr_init2(T60Q[k], PREC); mpfr_init2(QA[k], PREC);
        mpfr_init2(QB[k], PREC); mpfr_init2(QC[k], PREC);
    }
    for (int t = 0; t < 2 * MAXDEG + 2; t++)
        for (int k = 0; k < 4; k++) { mpfr_init2(PRE[t][k], PREC); mpfr_init2(SUF[t][k], PREC); }
    mpfr_inits2(PREC, CT1, CT2, CT3, (mpfr_ptr)0);
    for (int e = 0; e < MAXE; e++) mpfr_init2(BSAVE[e], PREC);
}

/* step quaternion for the corner angle: (cos(a/2), 0, 0, sin(a/2));
   called after ALPHA is set for the current solve */
static void set_step_quat(void)
{
    mpfr_div_ui(CT1, ALPHA, 2, MPFR_RNDN);
    mpfr_cos(T60Q[0], CT1, MPFR_RNDN);
    mpfr_set_ui(T60Q[1], 0, MPFR_RNDN); mpfr_set_ui(T60Q[2], 0, MPFR_RNDN);
    mpfr_sin(T60Q[3], CT1, MPFR_RNDN);
}

static void pair_word(int v, mpfr_t q[4])
{
    mpfr_set_ui(q[0], 1, MPFR_RNDN);
    for (int k = 1; k < 4; k++) mpfr_set_ui(q[k], 0, MPFR_RNDN);
    for (int t = 0; t < FLOWER_LEN[v]; t++) {
        q_mul(q, T60Q, q);
        int e = FLOWER_E[v][t];
        if (IS_FLAT[e]) continue;
        mpfr_set(QC[0], PC[e], MPFR_RNDN); mpfr_set(QC[1], PS[e], MPFR_RNDN);
        mpfr_set_ui(QC[2], 0, MPFR_RNDN);  mpfr_set_ui(QC[3], 0, MPFR_RNDN);
        q_mul(q, QC, q);
    }
}

static int cert_rows(void)   /* fills R2; returns row count */
{
    int r = 0;
    for (int k = 0; k < NACT; k++) {
        pair_word(ACT[k], QA);
        mpfr_sub_si(R2[r + 0], QA[0], SGN[ACT[k]], MPFR_RNDN);
        mpfr_set(R2[r + 1], QA[1], MPFR_RNDN);
        mpfr_set(R2[r + 2], QA[2], MPFR_RNDN);
        mpfr_set(R2[r + 3], QA[3], MPFR_RNDN);
        r += 4;
    }
    for (int f = 0; f < NFREE; f++) {
        int e = FREE_E[f];
        mpfr_mul(CT1, PC[e], PC[e], MPFR_RNDN);
        mpfr_mul(CT2, PS[e], PS[e], MPFR_RNDN);
        mpfr_add(CT1, CT1, CT2, MPFR_RNDN);
        mpfr_sub_ui(R2[r], CT1, 1, MPFR_RNDN);
        r++;
    }
    return r;
}

static void cert_jacobian(void)
{
    for (int e = 0; e < NE; e++)
        for (int s = 0; s < 2; s++)
            for (int c = 0; c < 4; c++)
                for (int k = 0; k < 2; k++) mpfr_set_ui(JB[e][s][c][k], 0, MPFR_RNDN);
    for (int kk = 0; kk < NACT; kk++) {
        int v = ACT[kk];
        /* factor list: T60 [pair] T60 [pair] ... */
        int nf = 0;
        static int slot_edge[MAXDEG]; int nslot = 0;
        static int factor_is_pair[2 * MAXDEG];
        static int factor_edge[2 * MAXDEG];
        for (int t = 0; t < FLOWER_LEN[v]; t++) {
            factor_is_pair[nf] = 0; factor_edge[nf] = -1; nf++;
            int e = FLOWER_E[v][t];
            if (!IS_FLAT[e]) {
                factor_is_pair[nf] = 1; factor_edge[nf] = e; nf++;
                slot_edge[nslot++] = e;
            }
        }
        (void)nslot;
        mpfr_set_ui(PRE[0][0], 1, MPFR_RNDN);
        for (int c = 1; c < 4; c++) mpfr_set_ui(PRE[0][c], 0, MPFR_RNDN);
        for (int t = 0; t < nf; t++) {
            if (factor_is_pair[t]) {
                int e = factor_edge[t];
                mpfr_set(QC[0], PC[e], MPFR_RNDN); mpfr_set(QC[1], PS[e], MPFR_RNDN);
                mpfr_set_ui(QC[2], 0, MPFR_RNDN);  mpfr_set_ui(QC[3], 0, MPFR_RNDN);
                q_mul(PRE[t], QC, PRE[t + 1]);
            } else q_mul(PRE[t], T60Q, PRE[t + 1]);
        }
        mpfr_set_ui(SUF[nf][0], 1, MPFR_RNDN);
        for (int c = 1; c < 4; c++) mpfr_set_ui(SUF[nf][c], 0, MPFR_RNDN);
        for (int t = nf - 1; t >= 0; t--) {
            if (factor_is_pair[t]) {
                int e = factor_edge[t];
                mpfr_set(QC[0], PC[e], MPFR_RNDN); mpfr_set(QC[1], PS[e], MPFR_RNDN);
                mpfr_set_ui(QC[2], 0, MPFR_RNDN);  mpfr_set_ui(QC[3], 0, MPFR_RNDN);
                q_mul(QC, SUF[t + 1], SUF[t]);
            } else q_mul(T60Q, SUF[t + 1], SUF[t]);
        }
        for (int t = 0; t < nf; t++) {
            if (!factor_is_pair[t]) continue;
            int e = factor_edge[t];
            int side = (EDGE_A[e] == v) ? 0 : 1;
            /* d/dc: pre * suf ; d/ds: pre * (0,1,0,0) * suf */
            q_mul(PRE[t], SUF[t + 1], QA);
            mpfr_set_ui(QC[0], 0, MPFR_RNDN); mpfr_set_ui(QC[1], 1, MPFR_RNDN);
            mpfr_set_ui(QC[2], 0, MPFR_RNDN); mpfr_set_ui(QC[3], 0, MPFR_RNDN);
            q_mul(PRE[t], QC, QB);
            q_mul(QB, SUF[t + 1], QB);
            for (int c = 0; c < 4; c++) {
                mpfr_add(JB[e][side][c][0], JB[e][side][c][0], QA[c], MPFR_RNDN);
                mpfr_add(JB[e][side][c][1], JB[e][side][c][1], QB[c], MPFR_RNDN);
            }
        }
    }
}

/* row structure helpers: active-vertex index per vertex id, fold cols */
static int ACTIDX[MAXV + 1];

static int row_nonzeros(int row, int *cols, double *valsd, mpfr_ptr *valsm)
{
    /* returns nonzero (column, value) entries of J row; values as mpfr
       pointers (valsm) and doubles (valsd). */
    int n = 0;
    if (row < 4 * NACT) {
        int k = row / 4, comp = row % 4, v = ACT[k];
        for (int t = 0; t < FLOWER_LEN[v]; t++) {
            int e = FLOWER_E[v][t];
            if (IS_FLAT[e]) continue;
            int side = (EDGE_A[e] == v) ? 0 : 1;
            int f = FIDX[e];
            cols[n] = 2 * f;     valsm[n] = JB[e][side][comp][0]; n++;
            cols[n] = 2 * f + 1; valsm[n] = JB[e][side][comp][1]; n++;
        }
    } else {
        int f = row - 4 * NACT;
        int e = FREE_E[f];
        mpfr_mul_ui(CT1, PC[e], 2, MPFR_RNDN);
        mpfr_mul_ui(CT2, PS[e], 2, MPFR_RNDN);
        cols[0] = 2 * f;     valsm[0] = CT1;
        cols[1] = 2 * f + 1; valsm[1] = CT2;
        n = 2;
    }
    for (int i = 0; i < n; i++) valsd[i] = mpfr_get_d(valsm[i], MPFR_RNDN);
    return n;
}

static int cert_chol(int n, mpfr_t M[MAXP][MAXP])
{
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < j; k++) {
            mpfr_mul(CT1, M[j][k], M[j][k], MPFR_RNDN);
            mpfr_sub(M[j][j], M[j][j], CT1, MPFR_RNDN);
        }
        if (mpfr_sgn(M[j][j]) <= 0) return -1;
        mpfr_sqrt(M[j][j], M[j][j], MPFR_RNDN);
        for (int i = j + 1; i < n; i++) {
            for (int k = 0; k < j; k++) {
                mpfr_mul(CT1, M[i][k], M[j][k], MPFR_RNDN);
                mpfr_sub(M[i][j], M[i][j], CT1, MPFR_RNDN);
            }
            mpfr_div(M[i][j], M[i][j], M[j][j], MPFR_RNDN);
        }
    }
    return 0;
}

static int certify(void)
{
    CERT_OK = 0; CERT_WHY = "";
    double U = ldexp(1.0, (int)(1 - PREC));
    /* actives, fold indices, signs, pairs from BEND */
    NACT = 0;
    for (int v = 1; v <= NV; v++) {
        ACTIDX[v] = -1;
        int any = 0;
        for (int t = 0; t < FLOWER_LEN[v]; t++)
            if (!IS_FLAT[FLOWER_E[v][t]]) { any = 1; break; }
        if (any) { ACTIDX[v] = NACT; ACT[NACT++] = v; }
    }
    for (int e = 0; e < NE; e++) FIDX[e] = -1;
    for (int f = 0; f < NFREE; f++) FIDX[FREE_E[f]] = f;
    int dmax = 0;
    for (int k = 0; k < NACT; k++)
        if (FLOWER_LEN[ACT[k]] > dmax) dmax = FLOWER_LEN[ACT[k]];
    double arity = 2.0 * dmax + 2.0;
    for (int f = 0; f < NFREE; f++) {
        int e = FREE_E[f];
        mpfr_div_ui(CT1, BEND[e], 2, MPFR_RNDN);
        mpfr_cos(PC[e], CT1, MPFR_RNDN);
        mpfr_sin(PS[e], CT1, MPFR_RNDN);
        mpfr_neg(PS[e], PS[e], MPFR_RNDN);      /* solver convention q_x(-b) */
    }
    for (int k = 0; k < NACT; k++) {
        pair_word(ACT[k], QA);
        SGN[ACT[k]] = (mpfr_sgn(QA[0]) >= 0) ? 1 : -1;
    }
    int ncol = 2 * NFREE;
    /* Newton polish: full-row normal equations, up to 4 iterations */
    for (int it = 0; it < 4; it++) {
        int nrows = cert_rows();
        double rn = 0;
        for (int i = 0; i < nrows; i++) {
            double x = fabs(mpfr_get_d(R2[i], MPFR_RNDN));
            if (x > rn) rn = x;
        }
        if (rn < 1e-33) break;
        cert_jacobian();
        for (int p = 0; p < ncol; p++) {
            mpfr_set_ui(G2[p], 0, MPFR_RNDN);
            for (int q = p; q < ncol; q++) mpfr_set_ui(A2[p][q], 0, MPFR_RNDN);
        }
        int cols[32]; double vd[32]; mpfr_ptr vm[32];
        for (int i = 0; i < nrows; i++) {
            int n = row_nonzeros(i, cols, vd, vm);
            for (int a = 0; a < n; a++) {
                mpfr_mul(CT3, vm[a], R2[i], MPFR_RNDN);
                mpfr_add(G2[cols[a]], G2[cols[a]], CT3, MPFR_RNDN);
                for (int b = a; b < n; b++) {
                    int p = cols[a] < cols[b] ? cols[a] : cols[b];
                    int q = cols[a] < cols[b] ? cols[b] : cols[a];
                    mpfr_mul(CT3, vm[a], vm[b], MPFR_RNDN);
                    mpfr_add(A2[p][q], A2[p][q], CT3, MPFR_RNDN);
                }
            }
        }
        for (int p = 0; p < ncol; p++)
            for (int q = 0; q < p; q++) mpfr_set(A2[p][q], A2[q][p], MPFR_RNDN);
        for (int p = 0; p < ncol; p++)
            for (int q = 0; q <= p; q++) mpfr_set(C2[p][q], A2[p][q], MPFR_RNDN);
        if (cert_chol(ncol, C2) != 0) { CERT_WHY = "newton chol failed"; return -1; }
        for (int i = 0; i < ncol; i++) {             /* solve C2 C2^T dx = g */
            mpfr_set(DX2[i], G2[i], MPFR_RNDN);
            for (int k = 0; k < i; k++) {
                mpfr_mul(CT1, C2[i][k], DX2[k], MPFR_RNDN);
                mpfr_sub(DX2[i], DX2[i], CT1, MPFR_RNDN);
            }
            mpfr_div(DX2[i], DX2[i], C2[i][i], MPFR_RNDN);
        }
        for (int i = ncol - 1; i >= 0; i--) {
            for (int k = i + 1; k < ncol; k++) {
                mpfr_mul(CT1, C2[k][i], DX2[k], MPFR_RNDN);
                mpfr_sub(DX2[i], DX2[i], CT1, MPFR_RNDN);
            }
            mpfr_div(DX2[i], DX2[i], C2[i][i], MPFR_RNDN);
        }
        for (int f = 0; f < NFREE; f++) {
            int e = FREE_E[f];
            mpfr_sub(PC[e], PC[e], DX2[2 * f], MPFR_RNDN);
            mpfr_sub(PS[e], PS[e], DX2[2 * f + 1], MPFR_RNDN);
        }
    }
    /* final residual + jacobian, row selection in doubles */
    int nrows = cert_rows();
    cert_jacobian();
    {
        int cols[32]; double vd[32]; mpfr_ptr vm[32];
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncol; j++) JD[i][j] = 0.0;
            int n = row_nonzeros(i, cols, vd, vm);
            for (int a = 0; a < n; a++) JD[i][cols[a]] = vd[a];
        }
    }
    {
        static int left[NRMAX]; int nleft = nrows;
        for (int i = 0; i < nrows; i++) { left[i] = i; INSEL[i] = 0; }
        for (int c = 0; c < ncol; c++) {
            int piv = -1; double pv = 1e-9;
            for (int li = 0; li < nleft; li++) {
                double x = fabs(JD[left[li]][c]);
                if (x > pv) { pv = x; piv = li; }
            }
            if (piv < 0) { CERT_WHY = "row selection failed"; return -1; }
            int prow = left[piv];
            SEL[c] = prow; INSEL[prow] = 1;
            left[piv] = left[--nleft];
            double inv = 1.0 / JD[prow][c];
            for (int li = 0; li < nleft; li++) {
                int r = left[li];
                double fac = JD[r][c] * inv;
                if (fac != 0.0)
                    for (int j = c; j < ncol; j++) JD[r][j] -= fac * JD[prow][j];
            }
        }
    }
    /* A = Js^T Js from selected rows (recomputed in mpfr) */
    for (int p = 0; p < ncol; p++)
        for (int q = p; q < ncol; q++) mpfr_set_ui(A2[p][q], 0, MPFR_RNDN);
    double amax = 0.0;
    {
        int cols[32]; double vd[32]; mpfr_ptr vm[32];
        for (int c = 0; c < ncol; c++) {
            int i = SEL[c];
            int n = row_nonzeros(i, cols, vd, vm);
            for (int a = 0; a < n; a++)
                for (int b = a; b < n; b++) {
                    int p = cols[a] < cols[b] ? cols[a] : cols[b];
                    int q = cols[a] < cols[b] ? cols[b] : cols[a];
                    mpfr_mul(CT3, vm[a], vm[b], MPFR_RNDN);
                    mpfr_add(A2[p][q], A2[p][q], CT3, MPFR_RNDN);
                }
        }
        for (int p = 0; p < ncol; p++)
            for (int q = 0; q < p; q++) mpfr_set(A2[p][q], A2[q][p], MPFR_RNDN);
        for (int p = 0; p < ncol; p++)
            for (int q = 0; q < ncol; q++) {
                double x = fabs(mpfr_get_d(A2[p][q], MPFR_RNDN));
                if (x > amax) amax = x;
            }
    }
    double EA = CSLOP * ncol * U * (amax + 1.0) * ncol;
    for (int p = 0; p < ncol; p++)
        for (int q = 0; q <= p; q++) mpfr_set(C2[p][q], A2[p][q], MPFR_RNDN);
    if (cert_chol(ncol, C2) != 0) { CERT_WHY = "cert chol failed"; return -1; }
    /* B = C^-1 (lower) */
    for (int i = 0; i < ncol; i++) {
        mpfr_ui_div(B2[i][i], 1, C2[i][i], MPFR_RNDN);
        for (int j = 0; j < i; j++) {
            mpfr_set_ui(CT2, 0, MPFR_RNDN);
            for (int k = j; k < i; k++) {
                mpfr_mul(CT1, C2[i][k], B2[k][j], MPFR_RNDN);
                mpfr_add(CT2, CT2, CT1, MPFR_RNDN);
            }
            mpfr_div(CT2, CT2, C2[i][i], MPFR_RNDN);
            mpfr_neg(B2[i][j], CT2, MPFR_RNDN);
        }
    }
    /* delta = ||B C - I||_F, bnorm = ||B||_F */
    mpfr_set_ui(CT2, 0, MPFR_RNDN);          /* dsum */
    for (int i = 0; i < ncol; i++)
        for (int j = 0; j <= i; j++) {
            mpfr_set_ui(CT3, 0, MPFR_RNDN);
            for (int k = j; k <= i; k++) {
                mpfr_mul(CT1, B2[i][k], C2[k][j], MPFR_RNDN);
                mpfr_add(CT3, CT3, CT1, MPFR_RNDN);
            }
            if (i == j) mpfr_sub_ui(CT3, CT3, 1, MPFR_RNDN);
            mpfr_mul(CT3, CT3, CT3, MPFR_RNDN);
            mpfr_add(CT2, CT2, CT3, MPFR_RNDN);
        }
    mpfr_sqrt(CT2, CT2, MPFR_RNDN);
    double delta = mpfr_get_d(CT2, MPFR_RNDU) + CSLOP * ncol * ncol * U;
    mpfr_set_ui(CT2, 0, MPFR_RNDN);
    for (int i = 0; i < ncol; i++)
        for (int j = 0; j <= i; j++) {
            mpfr_mul(CT1, B2[i][j], B2[i][j], MPFR_RNDN);
            mpfr_add(CT2, CT2, CT1, MPFR_RNDN);
        }
    mpfr_sqrt(CT2, CT2, MPFR_RNDN);
    double bnorm = mpfr_get_d(CT2, MPFR_RNDU);
    if (delta >= 1.0) { CERT_WHY = "witness delta>=1"; return -1; }
    double sigC = (1.0 - delta) / bnorm;
    /* rhoA = ||A - C C^T||_F + slops */
    mpfr_set_ui(CT2, 0, MPFR_RNDN);
    for (int i = 0; i < ncol; i++)
        for (int j = 0; j <= i; j++) {
            mpfr_set(CT3, A2[i][j], MPFR_RNDN);
            int kmax = (i < j ? i : j);
            for (int k = 0; k <= kmax; k++) {
                mpfr_mul(CT1, C2[i][k], C2[j][k], MPFR_RNDN);
                mpfr_sub(CT3, CT3, CT1, MPFR_RNDN);
            }
            mpfr_mul(CT3, CT3, CT3, MPFR_RNDN);
            if (i != j) mpfr_mul_ui(CT3, CT3, 2, MPFR_RNDN);
            mpfr_add(CT2, CT2, CT3, MPFR_RNDN);
        }
    mpfr_sqrt(CT2, CT2, MPFR_RNDN);
    double rhoA = mpfr_get_d(CT2, MPFR_RNDU) + CSLOP * ncol * ncol * U + EA;
    double margin = sigC * sigC - rhoA;
    if (margin <= 0.0) { CERT_WHY = "sigma margin<=0"; return -1; }
    double sig_low = sqrt(margin);
    /* selected residual norm */
    mpfr_set_ui(CT2, 0, MPFR_RNDN);
    for (int c = 0; c < ncol; c++) {
        mpfr_mul(CT1, R2[SEL[c]], R2[SEL[c]], MPFR_RNDN);
        mpfr_add(CT2, CT2, CT1, MPFR_RNDN);
    }
    mpfr_sqrt(CT2, CT2, MPFR_RNDN);
    double rnorm = mpfr_get_d(CT2, MPFR_RNDU) + CSLOP * ncol * arity * U;
    double beta = 1.0 / sig_low;
    double L = 4.0 * arity * arity;
    double eta = beta * rnorm;
    double h = beta * L * eta;
    if (!(h < 0.5)) { CERT_WHY = "kantorovich h>=1/2"; return -1; }
    /* stable form of (1 - sqrt(1-2h))/(beta L) for tiny h */
    double radius = 2.0 * eta / (1.0 + sqrt(1.0 - 2.0 * h));
    double drop_at = 0.0;
    for (int i = 0; i < nrows; i++) {
        if (INSEL[i]) continue;
        double x = fabs(mpfr_get_d(R2[i], MPFR_RNDN));
        if (x > drop_at) drop_at = x;
    }
    CERT_DROP_AT = drop_at;
    double drop = drop_at + CSLOP * arity * U + arity * radius;
    /* neoconvex window on the emitted bends, margin deg*2*radius */
    double tmin = 1e300, tmax = -1e300; int anyt = 0;
    for (int v = 1; v <= NV; v++) {
        if (ACTIDX[v] < 0) continue;
        double T = 0.0;
        for (int t = 0; t < FLOWER_LEN[v]; t++) {
            int e = FLOWER_E[v][t];
            if (!IS_FLAT[e]) T += mpfr_get_d(BEND[e], MPFR_RNDN);
        }
        double m = FLOWER_LEN[v] * 2.0 * radius + 1e-15;
        if (T - m < tmin) tmin = T - m;
        if (T + m > tmax) tmax = T + m;
        anyt = 1;
    }
    CERT_OK = 1; CERT_SIG = sig_low; CERT_H = h; CERT_RADIUS = radius;
    CERT_DROP = drop; CERT_ETA = eta;
    CERT_TMIN = anyt ? tmin : 0.0; CERT_TMAX = anyt ? tmax : 0.0;
    return 0;
}

/* the adopted record format (see FORMAT.md): halfturn units, certified-
   pair digits, integer tokens exact, proof internals as a # comment */
static int write_prove_record(FILE *fp, const char *name, const char *netcode)
{
    static mpfr_t PIC, BV, HT;
    static int init = 0;
    if (!init) { mpfr_inits2(PREC, PIC, BV, HT, (mpfr_ptr)0); init = 1; }
    mpfr_const_pi(PIC, MPFR_RNDN);
    fprintf(fp, "net %s\nv %d\ne %d\nunit halfturns\n", name, NV, NE);
    if (ALPHA_DEG != 60.0)
        fprintf(fp, "alpha %.17g degrees\n", ALPHA_DEG);
    if (PANCAKE) {
        fprintf(fp, "benderr 0\n");
        fprintf(fp, "faces %s\n", netcode);
        fprintf(fp, "# proof: degenerate pancake, all bends pinned exactly\n");
    } else {
        double r = CERT_OK ? CERT_RADIUS : -1.0;
        double benderr = CERT_OK ?
            (2.0 * r / (1.0 - r)) / 3.14159265358979 + 0.5e-36 : -1.0;
        fprintf(fp, "benderr %.1e\n", benderr);
        fprintf(fp, "faces %s\n", netcode);
        fprintf(fp, "# proof: kantorovich at %ld bits, jacobian sigma_min %.3e, "
                    "h %.3e, turning window [%.9e, %.9e], flats at %.0e, "
                    "solver euclid_lm_mp%s%s\n",
                (long)PREC, CERT_SIG, CERT_H, CERT_TMIN, CERT_TMAX,
                PROVE_TOL_USED,
                CERT_OK ? "" : ", CERTIFICATE FAILED: ",
                CERT_OK ? "" : CERT_WHY);
    }
    for (int a = 1; a <= NV; a++)
        for (int b = a + 1; b <= NV; b++) {
            int e = EDGE_IDX[a][b];
            if (e < 0) continue;
            if (IS_FLAT[e]) { fprintf(fp, "b %d %d 0\n", a, b); continue; }
            if (PANCAKE)    { fprintf(fp, "b %d %d %d\n", a, b, IS_PI[e]); continue; }
            /* certified-pair bend: pairs are (cos b/2, -sin b/2) */
            mpfr_atan2(BV, PS[e], PC[e], MPFR_RNDN);
            mpfr_mul_si(BV, BV, -2, MPFR_RNDN);
            mpfr_div(HT, BV, PIC, MPFR_RNDN);
            mpfr_fprintf(fp, "b %d %d %.36Rf\n", a, b, HT);
        }
    fprintf(fp, "end\n");
    return ferror(fp) ? -1 : 0;
}


/* ---------------- drivers ---------------- */

/* parse "a,b;c,d;..." into IS_FLAT; returns -1 on unknown pair */
static int parse_flats(const char *s)
{
    for (int e = 0; e < NE; e++) IS_FLAT[e] = 0;
    if (s) {
        const char *p = s;
        while (*p) {
            int a, b;
            if (sscanf(p, "%d,%d", &a, &b) != 2) return -1;
            int lo = a < b ? a : b, hi = a < b ? b : a;
            if (lo < 1 || hi > NV || EDGE_IDX[lo][hi] < 0) return -1;
            IS_FLAT[EDGE_IDX[lo][hi]] = 1;
            while (*p && *p != ';') p++;
            if (*p == ';') p++;
        }
    }
    NFREE = 0;
    for (int e = 0; e < NE; e++)
        if (!IS_FLAT[e]) FREE_E[NFREE++] = e;
    return 0;
}

/* seed file: lines "a b bend-decimal"; applies to free edges only */
static int load_seed(const char *path)
{
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;
    char buf[4096];
    while (fgets(buf, sizeof buf, fp)) {
        int a, b, off = 0;
        if (sscanf(buf, "%d %d %n", &a, &b, &off) < 2) continue;
        int lo = a < b ? a : b, hi = a < b ? b : a;
        if (lo < 1 || hi > NV || EDGE_IDX[lo][hi] < 0) { fclose(fp); return -1; }
        int e = EDGE_IDX[lo][hi];
        if (!IS_FLAT[e]) mpfr_set_str(BEND[e], buf + off, 10, MPFR_RNDN);
    }
    fclose(fp);
    return 0;
}

static int solve_one_netcode(const char *netcode, const char *name,
                             FILE *objout, char *errmsg, size_t errmsgsize)
{
    static double wish[MAXE];
    if (errmsgsize > 0) errmsg[0] = '\0';
    if (parse_netcode(netcode) <= 0) { snprintf(errmsg, errmsgsize, "parse failed"); return 1; }
    if (build_topology() < 0) { snprintf(errmsg, errmsgsize, "build_topology failed"); return 1; }
    mp_init_upto(NE, NV);
    cert_init_upto(2 * NE);
    if (parse_flats(FLATS_STR) < 0) { snprintf(errmsg, errmsgsize, "bad --flats"); return 1; }
    if (parse_dents(DENTS_STR) < 0) { snprintf(errmsg, errmsgsize, "bad --dents"); return 1; }
    if (wish_init_double(wish) < 0) { snprintf(errmsg, errmsgsize, "wish_init failed"); return 1; }
    for (int e = 0; e < NE; e++)
        mpfr_set_d(BEND[e], IS_FLAT[e] ? 0.0 : wish[e], MPFR_RNDN);
    if (SEED_PATH && load_seed(SEED_PATH) < 0) {
        snprintf(errmsg, errmsgsize, "bad --seed"); return 1;
    }
    mpfr_const_pi(ALPHA, MPFR_RNDN);
    mpfr_mul_d(ALPHA, ALPHA, ALPHA_DEG, MPFR_RNDN);
    mpfr_div_ui(ALPHA, ALPHA, 180, MPFR_RNDN);
    if (PROVE) set_step_quat();       /* cert buffers exist only under --prove */
    /* tol = 2^(-3*prec/4): ~1e-29 at 128 bits */
    mpfr_set_ui(TOL, 1, MPFR_RNDN);
    mpfr_div_2ui(TOL, TOL, (unsigned long)(3 * PREC / 4), MPFR_RNDN);
    if (lm_solve(200) < 0) { snprintf(errmsg, errmsgsize, "LM did not converge"); return 1; }
    if (PROVE) {
        /* pancake pattern: every bend within tol of 0 or +-gem/2 -> pin all
           exactly, no certificate needed (exact by construction) */
        PANCAKE = 1;
        double dpi = 3.14159265358979323846;
        for (int e = 0; e < NE; e++) {
            IS_PI[e] = 0;
            if (IS_FLAT[e]) continue;
            double b = mpfr_get_d(BEND[e], MPFR_RNDN);
            if (fabs(b) < PROVE_FLAT_TOL) continue;
            if (fabs(fabs(b) - dpi) < PROVE_FLAT_TOL) continue;
            PANCAKE = 0; break;
        }
        if (PANCAKE) {
            for (int e = 0; e < NE; e++) {
                double b = mpfr_get_d(BEND[e], MPFR_RNDN);
                if (!IS_FLAT[e] && fabs(b) < PROVE_FLAT_TOL) IS_FLAT[e] = 1;
                if (!IS_FLAT[e]) IS_PI[e] = (b > 0) ? 1 : -1;
            }
        } else {
            /* freeze ladder: classify at descending thresholds; accept only
               when the certificate passes AND the dropped closure rows are
               quiet (a wrongly frozen genuine fold shows up there -- the
               census tail of true folds reaches 4.6e-14, and correct
               solves measure drop_at <= 2.4e-26) */
            /* explicit --flats: trust the given flat set as the only rung
               (the ladder's threshold classification cannot separate exact
               zeros from genuine folds below 1e-8, e.g. the 00x00x lune
               nets whose fold tail reaches 1e-12) */
            static const double LADDER[4] = { 1e-8, 1e-11, 1e-14, 0.0 };
            for (int e = 0; e < NE; e++) mpfr_set(BSAVE[e], BEND[e], MPFR_RNDN);
            int accepted = 0;
            int nrungs = FLATS_STR ? 1 : 4;
            for (int li = 0; li < nrungs && !accepted; li++) {
                double tol = FLATS_STR ? -1 : LADDER[li];
                int newflats = 0;
                for (int e = 0; e < NE; e++) {
                    mpfr_set(BEND[e], BSAVE[e], MPFR_RNDN);
                    if (!FLATS_STR) IS_FLAT[e] = 0;
                }
                if (FLATS_STR) {
                    for (int e = 0; e < NE; e++)
                        if (IS_FLAT[e]) { mpfr_set_ui(BEND[e], 0, MPFR_RNDN); newflats++; }
                    PROVE_TOL_USED = -1;
                }
                for (int e = 0; e < NE; e++) {
                    if (tol > 0 && fabs(mpfr_get_d(BEND[e], MPFR_RNDN)) < tol) {
                        IS_FLAT[e] = 1;
                        mpfr_set_ui(BEND[e], 0, MPFR_RNDN);
                        newflats++;
                    }
                }
                NFREE = 0;
                for (int e = 0; e < NE; e++)
                    if (!IS_FLAT[e]) FREE_E[NFREE++] = e;
                if (newflats > 0 && NFREE > 0) {
                    if (lm_solve(60) < 0) {
                        if (getenv("LADDER_DEBUG"))
                            fprintf(stderr, "rung %.0e: %d flats, lm_solve FAILED\n", tol, newflats);
                        continue;      /* next rung */
                    }
                }
                int crc = certify();
                if (getenv("LADDER_DEBUG"))
                    fprintf(stderr, "rung %.0e: %d flats, certify rc=%d ok=%d drop_at=%.3e\n",
                            tol, newflats, crc, CERT_OK, CERT_DROP_AT);
                if (crc != 0) continue;
                if (CERT_OK && CERT_DROP_AT <= 1e-20) {
                    PROVE_TOL_USED = tol;
                    accepted = 1;
                }
            }
            if (!accepted) {
                snprintf(errmsg, errmsgsize, "certificate failed at every flat threshold");
                return 1;
            }
        }
    }
    if (PROVE) {
        if (write_prove_record(objout, name, netcode) < 0) {
            snprintf(errmsg, errmsgsize, "record write failed"); return 1;
        }
        return 0;
    }
    if (!BENDS_ONLY) {
        if (realize() < 0) { snprintf(errmsg, errmsgsize, "realize failed"); return 1; }
    }
    if (write_obj(objout) < 0) { snprintf(errmsg, errmsgsize, "OBJ write failed"); return 1; }
    return 0;
}

static void chomp(char *s)
{
    size_t n = strlen(s);
    while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r')) s[--n] = '\0';
}

int main(int argc, char **argv)
{
    int argi = 1;
    while (argi + 1 < argc) {
        if (strcmp(argv[argi], "--prec") == 0) {
            PREC = (mpfr_prec_t)atol(argv[argi + 1]);
            if (PREC < 24 || PREC > 8192) { fprintf(stderr, "bad prec\n"); return 2; }
            argi += 2;
        } else if (strcmp(argv[argi], "--flats") == 0) {
            FLATS_STR = argv[argi + 1]; argi += 2;
        } else if (strcmp(argv[argi], "--seed") == 0) {
            SEED_PATH = argv[argi + 1]; argi += 2;
        } else if (strcmp(argv[argi], "--sparse") == 0) {
            USE_SPARSE = 1; argi += 1;    /* default; kept for compatibility */
        } else if (strcmp(argv[argi], "--dense") == 0) {
            USE_SPARSE = 0; argi += 1;
        } else if (strcmp(argv[argi], "--nogate") == 0) {
            NOGATE = 1; argi += 1;
        } else if (strcmp(argv[argi], "--alpha") == 0) {
            ALPHA_DEG = atof(argv[argi + 1]); argi += 2;
            if (!(ALPHA_DEG > 0.0 && ALPHA_DEG <= 60.0)) {
                fprintf(stderr, "--alpha must be in (0, 60]\n");
                return 2;
            }
        } else if (strcmp(argv[argi], "--dents") == 0) {
            DENTS_STR = argv[argi + 1]; argi += 2;
        } else if (strcmp(argv[argi], "--dents-exact") == 0) {
            DENTS_EXACT = 1; argi += 1;
        } else if (strcmp(argv[argi], "--bends-only") == 0) {
            BENDS_ONLY = 1; argi += 1;
        } else if (strcmp(argv[argi], "--prove") == 0) {
            PROVE = 1; BENDS_ONLY = 1; argi += 1;
        } else if (strcmp(argv[argi], "--name") == 0) {
            NET_NAME = argv[argi + 1]; argi += 2;
        } else break;
    }
    while (argi < argc && (strcmp(argv[argi], "--bends-only") == 0 ||
                           strcmp(argv[argi], "--prove") == 0)) {
        if (strcmp(argv[argi], "--prove") == 0) { PROVE = 1; BENDS_ONLY = 1; }
        else BENDS_ONLY = 1;
        argi++;
    }
    mp_init_all();
    if (PROVE) cert_init_all();

    if (argi < argc && strcmp(argv[argi], "--batch") != 0) {
        char errmsg[512];
        int rc = solve_one_netcode(argv[argi], NET_NAME ? NET_NAME : "-", stdout, errmsg, sizeof(errmsg));
        if (rc != 0) { fprintf(stderr, "%s\n", errmsg[0] ? errmsg : "solver failed"); return 1; }
        return 0;
    }
    if (argi < argc && strcmp(argv[argi], "--batch") == 0) {
        /* batch: netcodes on stdin, one OBJ per line to stdout separated
           by "=== <index> <status>" headers (keep it simple: stream) */
        char *line = NULL; size_t cap = 0; ssize_t len; long idx = 0; int fails = 0;
        while ((len = getline(&line, &cap, stdin)) != -1) {
            (void)len; chomp(line);
            if (!line[0]) continue;
            char errmsg[512];
            char *nm = "-", *nc = line;
            char *sp = strchr(line, ' ');
            if (sp) { *sp = '\0'; nm = line; nc = sp + 1; }
            if (!PROVE) printf("=== %ld %s\n", idx, nc);
            int rc = solve_one_netcode(nc, nm, stdout, errmsg, sizeof(errmsg));
            if (!PROVE)
                printf("=== %ld %s %s\n", idx, rc == 0 ? "ok" : "fail", rc == 0 ? "" : errmsg);
            else if (rc != 0)
                printf("# failed %s %s\n", nc, errmsg);
            fails += (rc != 0);
            idx++;
        }
        free(line);
        return fails ? 1 : 0;
    }
    fprintf(stderr, "usage: %s [--prec BITS] [--flats \"a,b;c,d\"] [--seed FILE] [--dense] NETCODE | --batch < netcodes\n", argv[0]);
    return 2;
}

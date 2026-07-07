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

#define MAXV 128
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
    for (int i = 0; i < NV; i++) M[i][NV] = -2.0;
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
static mpfr_t A_BUF[MAXE][MAXE];  /* J^T J (+ lambda D on the copy) */
static mpfr_t M_BUF[MAXE][MAXE];
static mpfr_t G_BUF[MAXE], D_BUF[MAXE], DELTA[MAXE];
static mpfr_t VXYZ[MAXV + 1][3];
static mpfr_t ALPHA, T1, T2, T3, T4, NORM, NTRIAL, LAMBDA, TOL;
static Quat QTMP1, QTMP2, QSTEP, QD;
static Quat QS[MAXDEG], DQS[MAXDEG], PPRE[MAXDEG + 1], SSUF[MAXDEG + 1];
static mpfr_t W1[3], W2[3], W3[3], W4[3], W5[3], W6[3];

static void mp_init_all(void)
{
    for (int e = 0; e < MAXE; e++) {
        mpfr_init2(BEND[e], PREC); mpfr_init2(BEND_T[e], PREC);
        mpfr_init2(G_BUF[e], PREC); mpfr_init2(D_BUF[e], PREC); mpfr_init2(DELTA[e], PREC);
        for (int k = 0; k < 6; k++) mpfr_init2(JCOL[e][k], PREC);
    }
    for (int i = 0; i < 3 * MAXV; i++) { mpfr_init2(R_BUF[i], PREC); mpfr_init2(RT_BUF[i], PREC); }
    for (int p = 0; p < MAXE; p++)
        for (int q = 0; q < MAXE; q++) { mpfr_init2(A_BUF[p][q], PREC); mpfr_init2(M_BUF[p][q], PREC); }
    for (int v = 0; v <= MAXV; v++)
        for (int k = 0; k < 3; k++) mpfr_init2(VXYZ[v][k], PREC);
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

/* normal matrix A = J^T J via shared rows: columns p, q overlap only on
   rows of shared endpoint vertices. */
static void normal_matrix(void)
{
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

static int LM_ITERS = 0;

static int lm_solve(int maxiter)
{
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
        for (int pf = 0; pf < NFREE; pf++) mpfr_set(D_BUF[pf], A_BUF[pf][pf], MPFR_RNDN);

        int accepted = 0;
        for (int retry = 0; retry < 20; retry++) {
            if (chol_solve_neg() != 0) {
                mpfr_mul_ui(LAMBDA, LAMBDA, 10, MPFR_RNDN);
                if (mpfr_cmp_d(LAMBDA, 1e12) > 0) return -1;
                continue;
            }
            for (int e = 0; e < NE; e++) mpfr_set(BEND_T[e], BEND[e], MPFR_RNDN);
            for (int pf = 0; pf < NFREE; pf++)
                mpfr_add(BEND_T[FREE_E[pf]], BEND[FREE_E[pf]], DELTA[pf], MPFR_RNDN);
            int dent_v = has_dent(BEND_T);
            residual(BEND_T, RT_BUF);
            vec_norm2(NTRIAL, RT_BUF, rows);

            if (mpfr_cmp(NTRIAL, NORM) < 0 && dent_v == 0) {
                for (int e = 0; e < NE; e++) mpfr_set(BEND[e], BEND_T[e], MPFR_RNDN);
                for (int i = 0; i < rows; i++) mpfr_set(R_BUF[i], RT_BUF[i], MPFR_RNDN);
                mpfr_set(NORM, NTRIAL, MPFR_RNDN);
                mpfr_mul_d(LAMBDA, LAMBDA, 0.3, MPFR_RNDN);
                if (mpfr_cmp_d(LAMBDA, 1e-30) < 0) mpfr_set_d(LAMBDA, 1e-30, MPFR_RNDN);
                accepted = 1;
                break;
            }
            mpfr_mul_ui(LAMBDA, LAMBDA, 10, MPFR_RNDN);
            if (mpfr_cmp_d(LAMBDA, 1e12) > 0) return -1;
        }
        if (!accepted) return -1;
    }
    return (mpfr_cmp(NORM, TOL) <= 0) ? 0 : -1;
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
    for (int e = 0; e < NE; e++) {
        if (IS_FLAT[e])
            fprintf(fp, "# bend %d %d %d 0 flat\n", e, EDGE_A[e], EDGE_B[e]);
        else
            mpfr_fprintf(fp, "# bend %d %d %d %.40Re\n", e, EDGE_A[e], EDGE_B[e], BEND[e]);
    }
    for (int v = 1; v <= NV; v++)
        mpfr_fprintf(fp, "v %.40Rf %.40Rf %.40Rf\n",
                     VXYZ[v][0], VXYZ[v][1], VXYZ[v][2]);
    for (int fi = 0; fi < NF; fi++)
        fprintf(fp, "f %d %d %d\n", FACES[fi][0], FACES[fi][1], FACES[fi][2]);
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

static int solve_one_netcode(const char *netcode, FILE *objout,
                             char *errmsg, size_t errmsgsize)
{
    static double wish[MAXE];
    if (errmsgsize > 0) errmsg[0] = '\0';
    if (parse_netcode(netcode) <= 0) { snprintf(errmsg, errmsgsize, "parse failed"); return 1; }
    if (build_topology() < 0) { snprintf(errmsg, errmsgsize, "build_topology failed"); return 1; }
    if (parse_flats(FLATS_STR) < 0) { snprintf(errmsg, errmsgsize, "bad --flats"); return 1; }
    if (wish_init_double(wish) < 0) { snprintf(errmsg, errmsgsize, "wish_init failed"); return 1; }
    for (int e = 0; e < NE; e++)
        mpfr_set_d(BEND[e], IS_FLAT[e] ? 0.0 : wish[e], MPFR_RNDN);
    if (SEED_PATH && load_seed(SEED_PATH) < 0) {
        snprintf(errmsg, errmsgsize, "bad --seed"); return 1;
    }
    mpfr_const_pi(ALPHA, MPFR_RNDN);
    mpfr_div_ui(ALPHA, ALPHA, 3, MPFR_RNDN);
    /* tol = 2^(-3*prec/4): ~1e-29 at 128 bits */
    mpfr_set_ui(TOL, 1, MPFR_RNDN);
    mpfr_div_2ui(TOL, TOL, (unsigned long)(3 * PREC / 4), MPFR_RNDN);
    if (lm_solve(200) < 0) { snprintf(errmsg, errmsgsize, "LM did not converge"); return 1; }
    if (realize() < 0) { snprintf(errmsg, errmsgsize, "realize failed"); return 1; }
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
        } else break;
    }
    mp_init_all();

    if (argi < argc && strcmp(argv[argi], "--batch") != 0) {
        char errmsg[512];
        int rc = solve_one_netcode(argv[argi], stdout, errmsg, sizeof(errmsg));
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
            printf("=== %ld begin\n", idx);
            int rc = solve_one_netcode(line, stdout, errmsg, sizeof(errmsg));
            printf("=== %ld %s %s\n", idx, rc == 0 ? "ok" : "fail", rc == 0 ? "" : errmsg);
            fails += (rc != 0);
            idx++;
        }
        free(line);
        return fails ? 1 : 0;
    }
    fprintf(stderr, "usage: %s [--prec BITS] [--flats \"a,b;c,d\"] [--seed FILE] NETCODE | --batch < netcodes\n", argv[0]);
    return 2;
}

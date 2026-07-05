"""mpkernel.py -- take-1 (128-bit MPFR) merged quaternion solver + certificate.

Formulation (see the bendq program record): recognize ONLY flat. Each fold
edge carries quaternion pair (c,s) ~ (cos b/2, sin b/2) with unit row
c^2+s^2=1; vertex closure = quaternion word = eps_v * 1. Newton (float64
warm start optional; here direct 128-bit), then certificate: numerically
selected square row subsystem; Cholesky factor-witness sigma lower bound
(scheme of euclid_prover certify_sigma_lower); Kantorovich with explicit
multilinear Hessian constant.

CERTIFICATE SEMANTICS: existence + local uniqueness of an exact root of the
SELECTED SQUARE SUBSYSTEM within the reported radius, rounding errors
bounded by the standard 128-bit model with slop x64. Dropped closure rows
are bounded and reported; their discharge is the lune lemma + the monodromy
identity (the program's two structural obligations). Production version:
per-op directed rounding, principled row selection.
"""
import time
import gmpy2
from gmpy2 import mpfr

PREC = 128
gmpy2.get_context().precision = PREC
U = mpfr(2) ** (1 - PREC)
SLOP = 64

C30 = gmpy2.sqrt(mpfr(3)) / 2
ONE = mpfr(1); ZERO = mpfr(0)
T60 = (C30, ZERO, ZERO, mpfr(1)/2)

def qmul(a, b):
    return (a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3],
            a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
            a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
            a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0])

def vert_quat(star, flats, cs):
    q = (ONE, ZERO, ZERO, ZERO)
    for e in star:
        q = qmul(q, T60)
        if e in flats: continue
        c, s = cs[e]
        q = qmul(q, (c, s, ZERO, ZERO))
    return q

def residual_rows(act, stars, flats, folds, cs, signs):
    rows = []
    for vid in act:
        q = vert_quat(stars[vid], flats, cs)
        rows += [q[0] - signs[vid], q[1], q[2], q[3]]
    for e in folds:
        c, s = cs[e]
        rows.append(c*c + s*s - ONE)
    return rows

def jacobian(act, stars, flats, folds, cs):
    fidx = {e: i for i, e in enumerate(folds)}
    ncol = 2 * len(folds)
    rows = []
    for vid in act:
        factors = []; slots = []
        for e in stars[vid]:
            factors.append(T60)
            if e not in flats:
                c, s = cs[e]
                slots.append((e, len(factors)))
                factors.append((c, s, ZERO, ZERO))
        pre = [(ONE, ZERO, ZERO, ZERO)]
        for f in factors: pre.append(qmul(pre[-1], f))
        suf = [(ONE, ZERO, ZERO, ZERO)]
        for f in reversed(factors): suf.append(qmul(f, suf[-1]))
        suf.reverse()
        r4 = [[ZERO]*ncol for _ in range(4)]
        for e, k in slots:
            dc = qmul(pre[k], suf[k+1])
            ds = qmul(qmul(pre[k], (ZERO, ONE, ZERO, ZERO)), suf[k+1])
            j = fidx[e]
            for comp in range(4):
                r4[comp][2*j] = dc[comp]
                r4[comp][2*j+1] = ds[comp]
        rows += r4
    for e in folds:
        c, s = cs[e]
        r = [ZERO]*ncol
        j = fidx[e]
        r[2*j] = 2*c; r[2*j+1] = 2*s
        rows.append(r)
    return rows

def matTvec(J, v):
    n = len(J[0]); out = [ZERO]*n
    for i, row in enumerate(J):
        vi = v[i]
        if vi == 0: continue
        for j in range(n):
            if row[j] != 0: out[j] += row[j]*vi
    return out

def normal_matrix(J):
    n = len(J[0])
    A = [[ZERO]*n for _ in range(n)]
    for row in J:
        nz = [(j, row[j]) for j in range(n) if row[j] != 0]
        for a in range(len(nz)):
            ja, va = nz[a]
            for b in range(a, len(nz)):
                jb, vb = nz[b]
                A[ja][jb] += va*vb
    for i in range(n):
        for j in range(i): A[i][j] = A[j][i]
    return A

def chol(A):
    n = len(A)
    C = [[ZERO]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            s = A[i][j]
            for k in range(j): s -= C[i][k]*C[j][k]
            if i == j:
                if s <= 0: return None
                C[i][i] = gmpy2.sqrt(s)
            else:
                C[i][j] = s / C[j][j]
    return C

def tri_inv_lower(C):
    n = len(C)
    Bm = [[ZERO]*n for _ in range(n)]
    for i in range(n):
        Bm[i][i] = ONE / C[i][i]
        for j in range(i):
            s = ZERO
            for k in range(j, i): s += C[i][k]*Bm[k][j]
            Bm[i][j] = -s / C[i][i]
    return Bm

def chol_solve(C, b):
    n = len(C)
    y = [ZERO]*n
    for i in range(n):
        s = b[i]
        for k in range(i): s -= C[i][k]*y[k]
        y[i] = s / C[i][i]
    x = [ZERO]*n
    for i in range(n-1, -1, -1):
        s = y[i]
        for k in range(i+1, n): s -= C[k][i]*x[k]
        x[i] = s / C[i][i]
    return x

def fro(M):
    s = ZERO
    for row in M:
        for x in row: s += x*x
    return gmpy2.sqrt(s)

def select_rows(J):
    m, n = len(J), len(J[0])
    Fw = [[float(x) for x in row] for row in J]
    chosen = []; rows_left = list(range(m))
    for c in range(n):
        piv, pv = None, 1e-9
        for i in rows_left:
            if abs(Fw[i][c]) > pv: piv, pv = i, abs(Fw[i][c])
        if piv is None: return None
        chosen.append(piv); rows_left.remove(piv)
        inv = 1.0/Fw[piv][c]
        for i in rows_left:
            f = Fw[i][c]*inv
            if f != 0.0:
                ri = Fw[i]; rp = Fw[piv]
                for j in range(c, n): ri[j] -= f*rp[j]
    return chosen

def solve_and_certify(stars, flats, folds, seed_bends, iters=4):
    """seed_bends: eid -> float/str bend. Returns (cs, cert dict)."""
    t0 = time.time()
    act = [vid for vid in sorted(stars) if any(e not in flats for e in stars[vid])]
    cs = {}
    for e in folds:
        b = mpfr(str(seed_bends[e]))
        cs[e] = (gmpy2.cos(b/2), gmpy2.sin(b/2))     # seed only
    signs = {}
    for vid in act:
        q = vert_quat(stars[vid], flats, cs)
        signs[vid] = ONE if q[0] > 0 else -ONE
    for it in range(iters):
        r = residual_rows(act, stars, flats, folds, cs, signs)
        rn = max(abs(x) for x in r)
        if rn < mpfr(10) ** -33: break
        J = jacobian(act, stars, flats, folds, cs)
        A = normal_matrix(J)
        g = matTvec(J, r)
        C = chol(A)
        if C is None: return None, dict(ok=False, why='newton chol failed')
        dx = chol_solve(C, g)
        for i, e in enumerate(folds):
            c, s = cs[e]
            cs[e] = (c - dx[2*i], s - dx[2*i+1])
    F = len(folds); ncol = 2*F
    r_all = residual_rows(act, stars, flats, folds, cs, signs)
    J_all = jacobian(act, stars, flats, folds, cs)
    sel = select_rows(J_all)
    if sel is None: return cs, dict(ok=False, why='row selection failed')
    selset = set(sel)
    Js = [J_all[i] for i in sel]
    rs = [r_all[i] for i in sel]
    dmax = max(len(stars[v]) for v in act)
    arity = 2*dmax + 2
    A = normal_matrix(Js)
    EA = SLOP * ncol * float(U) * (max(max(abs(float(x)) for x in row) for row in A) + 1) * ncol
    C = chol(A)
    if C is None: return cs, dict(ok=False, why='cert chol failed')
    Bm = tri_inv_lower(C)
    n = ncol
    dsum = ZERO
    for i in range(n):
        for j in range(i+1):
            s = ZERO
            for k in range(j, i+1): s += Bm[i][k]*C[k][j]
            if i == j: s -= ONE
            dsum += s*s
    delta = gmpy2.sqrt(dsum) + SLOP*n*n*U
    bnorm = fro(Bm)
    if delta >= 1: return cs, dict(ok=False, why='witness delta>=1')
    sigC = (1 - delta)/bnorm
    rsum = ZERO
    for i in range(n):
        for j in range(i+1):
            s = A[i][j]
            for k in range(min(i, j)+1): s -= C[i][k]*C[j][k]
            rsum += s*s * (1 if i == j else 2)
    rhoA = gmpy2.sqrt(rsum) + SLOP*n*n*U + mpfr(EA)
    margin = sigC*sigC - rhoA
    if margin <= 0: return cs, dict(ok=False, why='sigma margin<=0')
    sig_low = gmpy2.sqrt(margin)
    rnorm = gmpy2.sqrt(sum(x*x for x in rs)) + SLOP*len(rs)*arity*U
    beta = 1/sig_low
    L = mpfr(4)*(arity*arity)
    eta = beta*rnorm
    h = beta*L*eta
    ok = h < mpfr(1)/2
    radius = (1 - gmpy2.sqrt(1 - 2*h))/(beta*L) if ok else mpfr('inf')
    drop_at = max((abs(r_all[i]) for i in range(len(r_all)) if i not in selset), default=ZERO)
    drop_bound = drop_at + SLOP*arity*U + mpfr(arity)*radius
    cert = dict(ok=bool(ok), F=F, nrows=len(r_all), sig_low=float(sig_low),
                eta=float(eta), h=float(h), radius=float(radius),
                drop_bound=float(drop_bound), secs=round(time.time()-t0, 2))
    return cs, cert

def bend_strings(cs, folds, digits=40):
    """display bends b = 2*atan2(s,c), high-precision decimal strings."""
    out = {}
    for e in folds:
        c, s = cs[e]
        b = 2*gmpy2.atan2(s, c)
        out[e] = gmpy2.digits(b, 10, digits)
    return out

def bend_str_format(b_digits):
    d, exp, _ = b_digits
    neg = d.startswith('-')
    if neg: d = d[1:]
    mant = d[0] + '.' + d[1:]
    return f"{'-' if neg else ''}{mant}e{exp-1}"

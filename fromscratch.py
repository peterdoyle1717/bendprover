"""fromscratch.py -- prototype: solve a net's bends from the net alone,
uniformly at 128 bits. No DP stage, no obj seeds, no dent gate.

    wish   : minimum-norm bend vector with vertex sums = gem (2*pi)
             (base-free: any vertex-share base projects to the same point)
    solve  : damped Gauss-Newton on the quaternion pair system (mpkernel
             residuals/Jacobian), Marquardt damping, no gate
    pop    : after convergence, if any vertex turn < 0, flip the signs of
             the (c,s) pairs on its incident fold edges and re-converge
             ("an ounce of cure")
    flats  : classify at FLAT_TOL on the converged bends, freeze, and hand
             to mpkernel.solve_and_certify for the frozen solve + take-1
             certificate

Driver (__main__): run over the census dumps, using each dump only for
the CLERS string (topology via the clers binary) and as the answer key.
"""
import math, sys, time
import gmpy2
from gmpy2 import mpfr
import mpkernel as K
import netio

FLAT_TOL = 1e-6      # post-MP-convergence flats sit ~1e-20; census min
                     # nonzero |bend| is 3.7e-4, population 3.8e-5
PI = math.pi

def wish_bends(eab, stars):
    """min-norm x with sum of x over each star = 1 rev; returns radians."""
    verts = sorted(stars); V = len(verts)
    vidx = {v: i for i, v in enumerate(verts)}
    edges = sorted(eab)
    # (B B^T) lam = -2c ; x = -B^T lam / 2   (base 0)
    M = [[0.0]*(V+1) for _ in range(V)]
    for v in verts: M[vidx[v]][vidx[v]] = float(len(stars[v]))
    for e, (a, b) in eab.items():
        M[vidx[a]][vidx[b]] += 1.0; M[vidx[b]][vidx[a]] += 1.0
    for i in range(V): M[i][V] = -2.0
    for c in range(V):
        p = max(range(c, V), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        for r in range(V):
            if r != c and M[r][c]:
                f = M[r][c]/M[c][c]
                for k in range(c, V+1): M[r][k] -= f*M[c][k]
    lam = [M[i][V]/M[i][i] for i in range(V)]
    return {e: -math.pi*(lam[vidx[a]] + lam[vidx[b]]) for e, (a, b) in eab.items()}

def rnorm(r):
    s = K.ZERO
    for x in r: s += x*x
    return gmpy2.sqrt(s)

def turns(stars, folds, cs):
    t = {}
    fold = set(folds)
    for v, star in stars.items():
        T = 0.0
        any_fold = False
        for e in star:
            if e in fold:
                c, s = cs[e]
                T += 2*math.atan2(float(s), float(c)); any_fold = True
        if any_fold: t[v] = T
    return t

def cur_signs(stars, cs):
    """per-vertex nearest target: sign of the current word's scalar part."""
    signs = {}
    for v in sorted(stars):
        q = K.vert_quat(stars[v], set(), cs)
        signs[v] = K.ONE if q[0] > 0 else -K.ONE
    return signs

def damped_gn(stars, folds, cs, maxit=60, lam0=1e-3):
    """Marquardt-damped GN on the full (vertex + unit) system. Returns
    (cs, resid, iters) -- no gate, no projection. Sign targets are
    re-derived from the current iterate each iteration (chase the
    nearest of +-1). Stops at 1e-12: the certify phase's undamped
    Newton finishes the job."""
    act = [v for v in sorted(stars)]
    flats = set()
    lam = mpfr(lam0)
    it = 0; slow = 0
    signs = cur_signs(stars, cs)
    r = K.residual_rows(act, stars, flats, folds, cs, signs)
    rn = rnorm(r)
    while it < maxit and rn > mpfr(10)**-12 and slow < 6:
        J = K.jacobian(act, stars, flats, folds, cs)
        A = K.normal_matrix(J)
        g = K.matTvec(J, r)
        accepted = False
        for _ in range(25):
            Ad = [row[:] for row in A]
            for i in range(len(Ad)): Ad[i][i] += lam*max(A[i][i], mpfr('1e-12'))
            C = K.chol(Ad)
            if C is None:
                lam *= 10; continue
            dx = K.chol_solve(C, g)
            cs2 = {}
            for i, e in enumerate(folds):
                c, s = cs[e]
                cs2[e] = (c - dx[2*i], s - dx[2*i+1])
            signs2 = cur_signs(stars, cs2)
            r2 = K.residual_rows(act, stars, flats, folds, cs2, signs2)
            rn2 = rnorm(r2)
            if rn2 < rn:
                slow = slow + 1 if rn2 > rn*mpfr('0.9') else 0
                cs, r, rn, signs = cs2, r2, rn2, signs2
                lam = max(lam*mpfr(0.3), mpfr('1e-14'))
                accepted = True
                break
            lam *= 10
            if lam > mpfr('1e15'): break
        if not accepted: break
        it += 1
    return cs, rn, it

def solve_net(eab, stars, max_pop_rounds=3, max_restarts=4):
    """from-scratch solve; returns dict with cs, bends, diagnostics.
    Wrong-basin cure: a converged root with pi edges on a mixed (not
    all-{0,pi}) configuration triggers a restart from a jittered wish."""
    import random
    t0 = time.time()
    folds = sorted(eab)
    w = wish_bends(eab, stars)
    rng = random.Random(len(eab)*1000003 + len(stars))
    pops = 0; rounds = 0; restarts = 0; rn = None
    for restarts in range(max_restarts+1):
        jit = 0.0 if restarts == 0 else 0.15
        cs = {e: (gmpy2.cos(mpfr(w[e]*(1+jit*(rng.random()-0.5)))/2),
                  gmpy2.sin(mpfr(w[e]*(1+jit*(rng.random()-0.5)))/2))
              for e in folds}
        good = False
        for rounds in range(1, max_pop_rounds+1):
            cs, rn, iters = damped_gn(stars, folds, cs)
            if rn > mpfr(10)**-10:
                break                   # no convergence: restart
            t = turns(stars, folds, cs)
            dented = [v for v, T in t.items() if T < -1e-9]
            if not dented:
                good = True             # converged, undented
                break
            if rounds == max_pop_rounds:
                break                   # dents persist: restart
            for v in dented:            # pop: mirror the star, re-converge
                for e in stars[v]:
                    c, s = cs[e]
                    cs[e] = (c, -s)
            pops += len(dented)
        if not good:
            continue
        bf = [abs(2*math.atan2(float(cs[e][1]), float(cs[e][0]))) for e in folds]
        npi = sum(1 for b in bf if abs(b - PI) < 1e-3)
        if npi == 0 or all(b < 1e-3 or abs(b - PI) < 1e-3 for b in bf):
            break                       # clean root, or genuine pancake
        # mixed pi/fold root on what should not be a pancake: wrong basin
    bends = {e: 2*gmpy2.atan2(cs[e][1], cs[e][0]) for e in folds}
    return dict(cs=cs, bends=bends, resid=rn, rounds=rounds, pops=pops,
                restarts=restarts, secs=round(time.time()-t0, 2))

def classify_and_certify(eab, stars, sol):
    """freeze flats found by the from-scratch solve, then run the
    production frozen solve + certificate (with the retry ladder)."""
    bends_f = {e: float(b) for e, b in sol['bends'].items()}
    pis = set(e for e, b in bends_f.items() if abs(abs(b) - PI) < 1e-3)
    if pis:
        # pancake (or unresolved wrong-basin root): never certified here
        flats = set(e for e, b in bends_f.items() if abs(b) < 1e-3)
        ok = (len(flats) + len(pis) == len(eab))
        return dict(ok=ok, pancake=True, flats=flats, pis=pis, folds=[],
                    why='' if ok else 'mixed pi/fold root after restarts')
    for tol in (1e-4, 1e-7):
        flats = set(e for e, b in bends_f.items() if abs(b) < tol)
        folds = sorted(e for e in eab if e not in flats)
        seed = {e: repr(bends_f[e]) for e in folds}
        cs, cert = K.solve_and_certify(stars, flats, folds, seed)
        cert['flats'] = flats; cert['pis'] = set(); cert['folds'] = folds
        cert['cs'] = cs; cert['flat_tol'] = tol
        if cert.get('ok'):
            return cert
    return cert

def compare_to_dump(eab, cert, sol, dump_bends, flat_tol=1e-4):
    """answer key: flats as sets; folds numerically with global mirror."""
    dflat = set(e for e in eab if abs(float(dump_bends[e])) < flat_tol)
    dpi = set(e for e in eab if abs(abs(float(dump_bends[e])) - PI) < flat_tol)
    ours = cert.get('flats', set())
    flats_match = (ours == dflat)
    if dpi:
        # pancake net: realizations form a degenerate family (tacos etc);
        # grade = recognized as pancake, all edges in {0, pi}
        return bool(cert.get('pancake') and cert.get('ok')), 0.0, flats_match
    if cert.get('pancake'):
        return False, float('nan'), flats_match   # pi root on non-pancake net
    if not cert.get('ok') or cert.get('cs') is None:
        return False, float('nan'), flats_match
    digs = K.bend_strings(cert['cs'], cert['folds'])
    ours_b = {e: float(K.bend_str_format(digs[e])) for e in cert['folds']}
    best = None
    for sgn in (1.0, -1.0):
        m = max(abs(sgn*ours_b[e] - float(dump_bends[e])) for e in cert['folds'])
        best = m if best is None else min(best, m)
    return cert.get('ok', False) and flats_match and best < 1e-9, best, flats_match

def main():
    import glob, os
    outp = sys.argv[1] if len(sys.argv) > 1 else "/tmp/fromscratch_322.tsv"
    dumps = sorted(glob.glob(
        "/Users/doyle/Dropbox/neo/_parking/bendq_sandbox/dumps/*.dump"))
    if len(sys.argv) > 2: dumps = dumps[:int(sys.argv[2])]
    npass = nfail = 0
    done = set()
    import os
    if os.path.exists(outp):           # resume: skip nets already graded
        for ln in open(outp):
            t = ln.split("\t")
            if t and t[0] != "clers": done.add(t[0])
    with open(outp, "a") as out:
        if not done:
            out.write("clers\tV\tstatus\trounds\tpops\titers_resid\tsecs\tmaxdiff\tflats_ok\trestarts\n")
        for p in dumps:
            clers, eab, dump_bends, stars = netio.read_dump(p)
            if clers in done: continue
            try:
                sol = solve_net(eab, stars)
                cert = classify_and_certify(eab, stars, sol)
                ok, maxdiff, flats_ok = compare_to_dump(eab, cert, sol, dump_bends)
            except Exception as ex:
                ok, maxdiff, flats_ok = False, float('nan'), False
                sol = dict(resid=-1, rounds=0, pops=0, restarts=0, secs=0)
                cert = dict(why=repr(ex)[:60])
            npass += ok; nfail += (not ok)
            out.write(f"{clers}\t{len(stars)}\t{'PASS' if ok else 'FAIL:'+str(cert.get('why','cmp'))[:40]}\t"
                      f"{sol['rounds']}\t{sol['pops']}\t{float(sol['resid']):.1e}\t"
                      f"{sol['secs']}\t{maxdiff:.2e}\t{flats_ok}\t{sol.get('restarts',0)}\n")
            out.flush()
            print(f"{'PASS' if ok else 'FAIL'}\t{clers[:36]}\tpops={sol['pops']} "
                  f"restarts={sol.get('restarts',0)} secs={sol['secs']} maxdiff={maxdiff:.1e}", flush=True)
    print(f"done: {npass} pass, {nfail} fail -> {outp}")

if __name__ == "__main__":
    main()

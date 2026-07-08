#!/usr/bin/env python3
"""embcheck.py -- certified embeddedness from a (bends, radius) record.

Input: certbends records ("=== idx NETCODE radius R" then lines
"a b bend"). The exact realization's vertices lie inside the ball
enclosures developed here (ball arithmetic over bend balls [b +- R]),
so every certificate below covers the exact polyhedron.

Tests, per the modular design:
  EDGE   shared-edge faces cross iff |bend| = gem/2; certify
         |b| + R < gem/2 ... i.e. hi(|b|) < pi.
  LINK   faces sharing exactly a vertex cross near it iff the vertex
         link polygon self-intersects on S^2. Certify simplicity:
         non-adjacent arc pairs disjoint (same-side fast path via the
         hemisphere lemma, else circle-intersection exclusion);
         adjacent arcs are covered by EDGE (fold-back needs |b| = pi).
  PAIR   vertex-disjoint simplex pairs (V+E+F): bounding-ball gap fast
         path, else candidate-axis separation, all in ball arithmetic.

LEMMA DEBT (documented, not yet written up): the reduction of shared-
vertex face crossings to LINK + vertex-disjoint incidences (the
plane-through-v argument, including the coplanar-sector case).

usage: embcheck.py RECORDS.certbends OUT.tsv
"""
import sys, math
from mpball import (Ball, vadd, vsub, vscale, vdot, vcross,
                    vnormpos, vunit)
import gmpy2
from gmpy2 import mpfr

PI = math.pi
COS_ARC = 0.5          # link arcs have length exactly pi/3

def parse_records(path):
    """FORMAT.md records: net/v/e/unit/benderr/faces/#/b.../end.
    Values are in halfturns; benderr is the per-bend bound in
    halfturns (0 = exact record)."""
    name = nc = None; benderr = 0.0; bends = {}
    for ln in open(path):
        t = ln.split()
        if not t or t[0] == "#": continue
        if t[0] == "net":
            name = t[1]; nc = None; benderr = 0.0; bends = {}
        elif t[0] == "benderr":
            benderr = float(t[1])
        elif t[0] == "faces":
            nc = t[1]
        elif t[0] == "b":
            bends[(int(t[1]), int(t[2]))] = t[3]
        elif t[0] == "end" and nc is not None:
            yield name, nc, benderr, bends
            name = nc = None; bends = {}

def bend_ball(valstr, benderr):
    """halfturn value string -> radian Ball. Integer tokens are exact;
    decimals carry benderr (halfturns)."""
    piB = Ball(gmpy2.const_pi(), 1e-55)
    if valstr in ("0",):
        return Ball(0, 0.0)
    if valstr in ("1", "-1"):
        return Ball(int(valstr)) * piB
    return Ball(mpfr(valstr), benderr) * piB

def develop_balls(faces, bend_of_pair, rad):
    """directed-edge BFS development in ball arithmetic; returns
    {vid: [Ball x3]} enclosing the exact realization's vertices."""
    apex = {}
    for a, b, c in faces:
        apex[(a, b)] = c; apex[(b, c)] = a; apex[(c, a)] = b
    a0, b0, c0 = faces[0]
    S3 = Ball(gmpy2.sqrt(mpfr(3)) / 2, 1e-56)
    pos = {a0: [Ball(0), Ball(0), Ball(0)],
           b0: [Ball(1), Ball(0), Ball(0)],
           c0: [Ball(0.5), S3, Ball(0)]}
    placed = {(a0, b0, c0)}
    from collections import deque
    dq = deque([(a0, b0, c0)])
    def fkey(f): return min(f, f[1:] + f[:1], f[2:] + f[:2])
    seen = {fkey((a0, b0, c0))}
    while dq:
        f0 = dq.popleft()
        for i in range(3):
            x, y = f0[i], f0[(i + 1) % 3]
            d = apex[(y, x)]
            f1 = (y, x, d)
            if fkey(f1) in seen: continue
            A, B, P = pos[x], pos[y], pos[f0[(i + 2) % 3]]
            # place_third (euclid_lm_mp form): additive radius growth,
            # no mirror-doubling of correlated errors
            M = vscale(Ball(0.5), vadd(A, B))
            e = vsub(B, A)
            w = vunit(vsub(P, M))
            n = vunit(vcross(w, e))
            th = bend_ball(bend_of_pair[tuple(sorted((x, y)))], rad)
            ct, st = th.cos(), th.sin()
            cperp = vadd(vscale(-ct, w), vscale(st, n))
            D = vadd(M, vscale(Ball(gmpy2.sqrt(mpfr(3)) / 2, 1e-56), cperp))
            pos[d] = D
            seen.add(fkey(f1)); dq.append(f1)
    return pos

def run_normal(run):
    """great-circle normal of a merged arc (consecutive dirs 60 deg apart,
    so the first sub-arc's cross product is bounded away from zero)."""
    return vcross(run[0], run[1])

def run_one_side(n, run):
    """whole merged arc strictly on one side of the circle with normal n:
    per sub-arc hemisphere lemma (sub-arcs have length pi/3 < pi)."""
    los = [vdot(n, d).lo() for d in run]
    his = [vdot(n, d).hi() for d in run]
    return all(l > 0 for l in los) or all(h < 0 for h in his)

def out_of_run(x, run):
    """certify x is outside every sub-arc of the merged arc."""
    for i in range(len(run) - 1):
        if not (vdot(x, run[i]).hi() < COS_ARC or
                vdot(x, run[i + 1]).hi() < COS_ARC):
            return False
    return True

def merged_pair_disjoint(P, Q):
    """certify two merged arcs (lists of dirs) disjoint."""
    nP, nQ = run_normal(P), run_normal(Q)
    if run_one_side(nP, Q) or run_one_side(nQ, P):
        return True
    # distinct circles meet at +-x only; exclude x from one arc
    xr = vcross(nP, nQ)
    if not vdot(xr, xr).lo() > 0:
        return False             # cannot certify circles distinct
    x = vunit(xr)
    for xx in (x, [-c for c in x]):
        if not (out_of_run(xx, P) or out_of_run(xx, Q)):
            return False
    return True

def check_net(nc, rad, bends):
    """rad and bend values in halfturns (FORMAT.md)."""
    faces = [tuple(int(x) for x in f.split(",")) for f in nc.split(";")]
    verts = sorted(set(v for f in faces for v in f))
    # exactly-flat records are pancakes: certified by exactness, not here
    if all(s in ("0", "1", "-1") for s in bends.values()):
        return "PANCAKE", "all bends integer halfturns (flat pancake, not embedded)"
    # EDGE test: shared-edge faces cross iff |bend| = 1 halfturn
    for (a, b), s in bends.items():
        if s in ("1", "-1") or abs(float(s)) + rad >= 1.0 - 1e-14:
            return "FAIL", f"edge ({a},{b}) bend at gem/2"
    pos = develop_balls(faces, bends, rad)
    # stars in cyclic order
    apex = {}
    for a, b, c in faces:
        apex[(a, b)] = c; apex[(b, c)] = a; apex[(c, a)] = b
    for v in verts:
        nb0 = next(w for (x, w) in apex if x == v)
        cyc = [nb0]
        while True:
            w = apex[(v, cyc[-1])]
            if w == nb0: break
            cyc.append(w)
        d = len(cyc)
        # exact flats from the record: corner at cyc[i] is flat iff the
        # edge (v, cyc[i]) has bend identically 0 -> its two faces are
        # exactly coplanar and the arcs merge structurally
        flatc = [bends[tuple(sorted((v, cyc[i])))] == "0" for i in range(d)]
        if all(flatc):
            continue             # fully flat vertex: star is a flat disk
        dirs = [vunit(vsub(pos[w], pos[v])) for w in cyc]
        start = next(i for i in range(d) if not flatc[i])
        runs = []; cur = [dirs[start]]
        for k in range(1, d + 1):
            i = (start + k) % d
            cur.append(dirs[i])
            if not flatc[i]:
                runs.append(cur); cur = [dirs[i]]
        m = len(runs)
        if m == 1:
            return "FAIL", f"single-fold link at vertex {v} (cannot certify)"
        if m == 2:
            continue             # lune: two arcs joined at both corners;
                                 # covered by the EDGE test (folds != gem/2)
        for i in range(m):
            for j in range(i + 1, m):
                if (j - i) % m in (1, m - 1): continue   # adjacent merged arcs
                if not merged_pair_disjoint(runs[i], runs[j]):
                    return "FAIL", f"link not certified simple at vertex {v} (merged arcs {i},{j})"
    # PAIR test: simplices, ball gaps, candidate axes
    simps = ([(v,) for v in verts]
             + [tuple(sorted(p)) for p in bends]
             + [tuple(f) for f in faces])
    cents = {}; rads = {}
    for s in simps:
        c = pos[s[0]]
        for w in s[1:]:
            c = vadd(c, pos[w])
        c = vscale(Ball(1.0 / len(s)), c)
        R = 0.0
        for w in s:
            dv = vsub(pos[w], c)
            R = max(R, math.sqrt(max(vdot(dv, dv).hi(), 0.0)))
        cents[s] = c; rads[s] = R
    n = len(simps)
    for i in range(n):
        si = simps[i]
        for j in range(i + 1, n):
            sj = simps[j]
            if set(si) & set(sj): continue
            dc = vsub(cents[si], cents[sj])
            gap2 = vdot(dc, dc)
            need = rads[si] + rads[sj]
            if gap2.lo() > need * need:
                continue                       # bounding-ball fast path
            ok = False
            cands = [dc]
            for a in si:
                for b in sj:
                    cands.append(vsub(pos[a], pos[b]))
            for nvec in cands:
                nn = vdot(nvec, nvec)
                if not nn.lo() > 0: continue
                lo_i = min(vdot(nvec, pos[a]).lo() for a in si)
                hi_j = max(vdot(nvec, pos[b]).hi() for b in sj)
                if lo_i > hi_j: ok = True; break
                hi_i = max(vdot(nvec, pos[a]).hi() for a in si)
                lo_j = min(vdot(nvec, pos[b]).lo() for b in sj)
                if lo_j > hi_i: ok = True; break
            if not ok:
                return "FAIL", f"separation not certified for {si} vs {sj}"
    return "PASS", ""

def main():
    recs, out_tsv = sys.argv[1], sys.argv[2]
    np_ = nf_ = 0
    with open(out_tsv, "w") as out:
        out.write("idx\tV\tembedded\twhy\n")
        for idx, nc, rad, bends in parse_records(recs):
            try:
                verdict, why = check_net(nc, rad, bends)
            except Exception as ex:
                verdict, why = "FAIL", repr(ex)[:70]
            V = max(x for f in nc.split(";") for x in
                    (int(y) for y in f.split(",")))
            np_ += (verdict != "FAIL"); nf_ += (verdict == "FAIL")
            out.write(f"{idx}\t{V}\t{verdict}\t{why}\n")
    print(f"done: {np_} embedded-or-pancake, {nf_} FAIL")

if __name__ == "__main__":
    main()

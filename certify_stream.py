#!/usr/bin/env python3
"""certify_stream.py -- certify existence from a euclid_lm_mp bends stream.

Reads the --bends-only --batch stream (records "=== idx NETCODE", "# bend
e a b VALUE" lines, "=== idx ok"). Per net:

  classify flats     |b| < 1e-8 on the 128-bit prelim bends (measured gap:
                     flats <= ~1e-15, genuine folds >= ~3.8e-5)
  frozen solve+cert  mpkernel.solve_and_certify: Newton on the frozen
                     quaternion system seeded by the prelim bends, then the
                     take-1 certificate (Cholesky witness + Kantorovich)
  neoconvex gate     certified root must have 0 < T_v < gem at every
                     non-flat vertex, with margin deg*radius; flat vertices
                     are exactly flat by construction

Emits one TSV line per net; for nets with flats, appends the certified
fold bends (exact flats) to a .certbends stream.

usage: certify_stream.py CHUNK.bends OUT.tsv [OUT.certbends] [--vmax N]
"""
import sys, math
import mpkernel as K
import netio

FLAT_TOL = 1e-8
GEM = 2 * math.pi

def parse_records(path):
    """yields (idx, netcode, {(a,b): bendstr}) for ok-terminated records."""
    idx = None; nc = None; bends = {}
    for ln in open(path):
        if ln.startswith("=== "):
            t = ln.split()
            if len(t) >= 3 and t[2] == "ok" and nc is not None:
                yield idx, nc, bends
                idx = nc = None; bends = {}
            elif len(t) == 3 and t[2] not in ("ok", "fail"):
                idx, nc, bends = t[1], t[2], {}
            elif len(t) >= 3 and t[2] == "fail":
                idx = nc = None; bends = {}
        elif nc is not None and ln.startswith("# bend"):
            t = ln.split()
            if len(t) >= 6:
                bends[(int(t[3]), int(t[4]))] = t[5]

def certify_one(nc, bendmap):
    faces = [tuple(int(x) for x in f.split(",")) for f in nc.split(";")]
    eab, stars = netio.stars_from_faces(faces)
    pair2e = {tuple(sorted(p)): e for e, p in eab.items()}
    seed = {}
    for (a, b), s in bendmap.items():
        seed[pair2e[tuple(sorted((a, b)))]] = s
    if len(seed) != len(eab):
        return dict(ok=False, why="bend/edge mismatch"), None, None
    flats = set(e for e, s in seed.items() if abs(float(s)) < FLAT_TOL)
    folds = sorted(e for e in eab if e not in flats)
    cs, cert = K.solve_and_certify(stars, flats, folds, seed)
    cert["V"] = len(stars); cert["nflat"] = len(flats)
    if not cert.get("ok"):
        return cert, None, None
    # certified fold bends and the neoconvex window gate
    digs = K.bend_strings(cs, folds)
    fout = {e: K.bend_str_format(digs[e]) for e in folds}
    radius = cert.get("radius", 0.0)
    tmin = None; tmax = None
    for v, star in stars.items():
        if all(e in flats for e in star):
            continue                     # exactly flat by construction
        T = sum(0.0 if e in flats else float(fout[e]) for e in star)
        margin = len(star) * radius
        lo, hi = T - margin, T + margin
        tmin = lo if tmin is None else min(tmin, lo)
        tmax = hi if tmax is None else max(tmax, hi)
    cert["tmin"] = tmin; cert["tmax"] = tmax
    # modular verdicts: ok = existence certificate; neoconvex = the
    # independent window check on the certified bends + radius
    cert["neoconvex"] = bool(tmin is not None and tmin > 0 and tmax < GEM)
    return cert, fout, (eab, flats)

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    vmax = None
    for i, a in enumerate(sys.argv):
        if a == "--vmax": vmax = int(sys.argv[i + 1])
    chunk, out_tsv = args[0], args[1]
    out_cb = open(args[2], "a") if len(args) > 2 else None
    npass = nfail = nskip = 0
    with open(out_tsv, "w") as out:
        out.write("idx\tV\tok\tneoconvex\tnflat\tsig_low\th\tradius\tdrop_bound\ttmin\ttmax\twhy\n")
        for idx, nc, bendmap in parse_records(chunk):
            try:
                V = max(x for f in nc.split(";") for x in
                        (int(y) for y in f.split(",")))
                if vmax is not None and V > vmax:
                    nskip += 1; continue
                cert, fout, extra = certify_one(nc, bendmap)
            except Exception as ex:
                cert = dict(ok=False, why=repr(ex)[:60], V=-1, nflat=-1)
                fout = None
            ok = bool(cert.get("ok"))
            npass += ok; nfail += (not ok)
            out.write(f"{idx}\t{cert.get('V')}\t{'PASS' if ok else 'FAIL'}\t"
                      f"{'Y' if cert.get('neoconvex') else 'N'}\t"
                      f"{cert.get('nflat')}\t{cert.get('sig_low','')}\t"
                      f"{cert.get('h','')}\t{cert.get('radius','')}\t"
                      f"{cert.get('drop_bound','')}\t{cert.get('tmin','')}\t"
                      f"{cert.get('tmax','')}\t{cert.get('why','')}\n")
            if ok and fout is not None and out_cb:
                # THE deliverable: bends (flats identically 0) + error bound
                eab, flats = extra
                out_cb.write(f"=== {idx} {nc} radius {cert.get('radius')}\n")
                for e in sorted(eab):
                    a, b = eab[e]
                    out_cb.write(f"{a} {b} {'0' if e in flats else fout[e]}\n")
    if out_cb: out_cb.close()
    print(f"done {chunk}: {npass} pass, {nfail} fail, {nskip} skipped")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""bendprover -- take-1 solver + certificate for unit-equilateral nets.

Per net: classify flats from a seed (flat gap in the v<=50 census: 3.7e-4,
threshold default 1e-4), pancake-detect (any bend near tau/2 -> degenerate
flat object, passthrough), solve the merged quaternion system at 128 bits,
certify (sigma factor witness + Kantorovich), and emit:

    OUT/v<NN><CLERS>.bends   exact 0 for flats, 40-digit folds, cert summary
    OUT/v<NN><CLERS>.cert    full certificate record
    OUT/v<NN><CLERS>.obj     developed coordinates (embeddedness input)

Seeds: --dump FILE (bend_sweep dump; CLERS read from its CASE line), or
--obj FILE (coords; bends derived, global mirror sign auto-fixed by closure
residual -- pass --clers to name the net). Batch: --batch LIST OUT where
LIST lines are "<seedpath>" or "<CLERS><TAB><seedpath>".

usage:
  bendprover.py [--clers CLERS] (--dump F | --obj F) OUTDIR
  bendprover.py --batch LIST OUTDIR
"""
import sys, os, math, json
import netio, mpkernel

FLAT_TOL = 1e-4
PI = math.pi

def classify(bends):
    flats = set(); pis = set(); folds = []
    for e, b in bends.items():
        fb = float(b)
        if abs(fb) < FLAT_TOL: flats.add(e)
        elif abs(abs(fb) - PI) < FLAT_TOL: pis.add(e)
        else: folds.append(e)
    return flats, pis, sorted(folds)

def closure_ok(stars, flats, bends):
    """max |H - Id| over active vertices at the seed (float, sign check)."""
    import math
    def rz(t):
        c, s = math.cos(t), math.sin(t); return [[c, -s, 0], [s, c, 0], [0, 0, 1]]
    def rx(t):
        c, s = math.cos(t), math.sin(t); return [[1, 0, 0], [0, c, -s], [0, s, c]]
    def mul(a, b): return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    worst = 0.0
    for vid, star in stars.items():
        if all(e in flats for e in star): continue
        H = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        for e in star:
            H = mul(H, rz(PI/3)); H = mul(H, rx(float(bends[e])))
        worst = max(worst, max(abs(H[i][j] - (1 if i == j else 0)) for i in range(3) for j in range(3)))
    return worst

def turning_report(stars, flats, bends_f):
    tmin = None; nflatv = 0
    for vid, star in stars.items():
        if all(e in flats for e in star): nflatv += 1; continue
        T = sum(bends_f[e] for e in star)
        tmin = T if tmin is None else min(tmin, T)
    return dict(min_nonflat=(round(tmin, 9) if tmin is not None else None), nflatv=nflatv)

def get_faces(clers, eab, clers_bin):
    """decode faces and check the vertex numbering agrees with the seed's
    edge pairs (same decoder upstream); return None if not available."""
    if not clers: return None
    faces = netio.faces_from_clers(clers, clers_bin)
    pairs = set(tuple(sorted((f[i], f[(i+1) % 3]))) for f in faces for i in range(3))
    if pairs != set(tuple(sorted(p)) for p in eab.values()):
        raise RuntimeError("clers faces disagree with seed edge pairs")
    return faces

def process(clers, seedpath, outdir, clers_bin="clers"):
    if seedpath.endswith(".dump"):
        dclers, eab, bends_s, stars = netio.read_dump(seedpath)
        clers = clers or dclers
    else:
        eab, bends_s, stars = netio.bends_from_obj(seedpath)
    faces = get_faces(clers, eab, clers_bin)
    bends_f = {e: float(b) for e, b in bends_s.items()}
    flats, pis, folds = classify(bends_s)
    # global mirror fix for obj seeds
    w = closure_ok(stars, flats, bends_f)
    if w > 1e-6:
        bends_f2 = {e: -b for e, b in bends_f.items()}
        if closure_ok(stars, flats, bends_f2) < w:
            bends_f = bends_f2
            bends_s = {e: repr(b) for e, b in bends_f.items()}
    nv = len(stars)
    stem = os.path.join(outdir, f"v{nv}{clers}" if clers else f"v{nv}seed_{os.path.basename(seedpath)}")
    if pis:
        # pancake: degenerate flat object (every edge 0 or tau/2); nothing to solve
        if folds:
            cert = dict(ok=False, pancake=True, why='pi edge outside a fully degenerate pancake')
        else:
            cert = dict(ok=True, pancake=True, F=0,
                        why='pancake: every edge 0 or tau/2; degenerate flat object, special-cased')
        tr = turning_report(stars, flats, bends_f)
        cert['turning'] = tr
        netio.write_bends(stem + ".bends", clers or "?", eab, flats,
                          {e: ("tau/2" if e in pis else bends_s[e])
                           for e in eab if e not in flats}, cert, tr)
        if faces:
            coords, worst, sign = netio.develop_best(faces, eab, bends_f)
            cert['develop_spread'] = worst
            netio.write_obj(stem + ".obj", faces, coords)
        with open(stem + ".cert", "w") as f: json.dump(cert, f, indent=1)
        return cert
    cs, cert = mpkernel.solve_and_certify(stars, flats, folds, bends_s)
    fold_out = {}
    if cs is not None:
        digs = mpkernel.bend_strings(cs, folds)
        fold_out = {e: mpkernel.bend_str_format(digs[e]) for e in folds}
        bends_exact = {e: 0.0 for e in eab}
        for e in folds: bends_exact[e] = float(fold_out[e])
    else:
        bends_exact = {e: (0.0 if e in flats else bends_f[e]) for e in eab}
    tr = turning_report(stars, flats, bends_exact)
    cert['turning'] = tr
    netio.write_bends(stem + ".bends", clers or "?", eab, flats, fold_out or bends_s, cert, tr)
    if faces:
        coords, worst, sign = netio.develop_best(faces, eab, bends_exact)
        cert['develop_spread'] = worst
        netio.write_obj(stem + ".obj", faces, coords)
    with open(stem + ".cert", "w") as f: json.dump(cert, f, indent=1)
    return cert

def main():
    args = sys.argv[1:]
    clers_bin = os.environ.get("CLERS_BIN", "clers")
    if args and args[0] == "--batch":
        listfile, outdir = args[1], args[2]
        os.makedirs(outdir, exist_ok=True)
        npass = nfail = 0
        for ln in open(listfile):
            t = [x.strip() for x in ln.rstrip("\n").split("\t")]
            if not t or not t[0] or t[0].startswith("#"): continue
            # "seedpath" (CLERS from dump CASE line) or "clers<TAB>seedpath"
            clers, seed = (None, t[0]) if len(t) == 1 else (t[0], t[1])
            try:
                cert = process(clers, seed, outdir, clers_bin)
                ok = cert.get('ok')
            except Exception as ex:
                ok = False; cert = dict(ok=False, why=str(ex)[:80])
            npass += bool(ok); nfail += (not ok)
            tag = clers or os.path.basename(seed)
            print(f"{'PASS' if ok else 'FAIL'}\t{tag[:44]}\t"
                  f"{cert.get('sig_low', cert.get('why', ''))}", flush=True)
        print(f"done: {npass} pass, {nfail} fail")
    else:
        clers = None; dump = None; obj = None; outdir = None
        it = iter(args)
        for a in it:
            if a == "--clers": clers = next(it)
            elif a == "--dump": dump = next(it)
            elif a == "--obj": obj = next(it)
            else: outdir = a
        if not outdir or not (dump or obj):
            print(__doc__); sys.exit(1)
        os.makedirs(outdir, exist_ok=True)
        cert = process(clers, dump or obj, outdir, clers_bin)
        print(json.dumps(cert, indent=1))

if __name__ == "__main__":
    main()

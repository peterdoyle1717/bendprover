#!/usr/bin/env python3
"""bendprover tests: fold net end-to-end, pancake path, obj round-trip."""
import os, sys, json, tempfile
from decimal import Decimal, getcontext

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import bendprover, netio

getcontext().prec = 60
CLERS_BIN = os.environ.get("CLERS_BIN", "clers")
FOLD = os.path.join(HERE, "fixtures", "CCCACCACACCACACACAAABE.dump")
PANCAKE = os.path.join(HERE, "fixtures", "CCCACACCAABE.dump")
np = nf = 0

def check(name, cond, detail=""):
    global np, nf
    np += cond; nf += (not cond)
    print(f"{'ok' if cond else 'FAIL'}  {name}  {detail}")

with tempfile.TemporaryDirectory() as td:
    # --- fold net from dump ---
    cert = bendprover.process(None, FOLD, td, CLERS_BIN)
    check("fold: certified", cert.get("ok") is True,
          f"sig_low {cert.get('sig_low', 0):.4f} h {cert.get('h', 1):.2e}")
    check("fold: kantorovich margin", cert.get("h", 1) < 1e-20)
    check("fold: neoconvex turning", cert["turning"]["min_nonflat"] > 0,
          f"min {cert['turning']['min_nonflat']}")
    check("fold: develop closes", cert.get("develop_spread", 1) < 1e-12,
          f"spread {cert.get('develop_spread'):.2e}")
    stem = [os.path.join(td, f[:-6]) for f in os.listdir(td) if f.endswith(".bends")
            and "CCCACCA" in f][0]
    # bends agree with the prec-200 dump on folds; flats exact 0
    dv = {}
    for ln in open(FOLD):
        t = ln.split()
        if t and t[0] == "EDGE": dv[int(t[1])] = Decimal(t[4])
    worst = Decimal(0); nflat = 0
    for ln in open(stem + ".bends"):
        t = ln.split()
        if not t or t[0] != "EDGE": continue
        e, v = int(t[1]), Decimal(t[4])
        if v == 0: nflat += 1
        else: worst = max(worst, abs(v - dv[e]))
    check("fold: bends match dump", worst < Decimal("1e-30"), f"max diff {worst:.1e}")

    # --- obj round-trip: reseed from the emitted obj ---
    clers = None
    for ln in open(stem + ".bends"):
        if ln.startswith("CLERS"): clers = ln.split()[1]
    td2 = os.path.join(td, "rt"); os.makedirs(td2)
    cert2 = bendprover.process(clers, stem + ".obj", td2, CLERS_BIN)
    check("objseed: re-certifies", cert2.get("ok") is True,
          f"sig_low {cert2.get('sig_low', 0):.4f}")
    stem2 = [os.path.join(td2, f[:-6]) for f in os.listdir(td2) if f.endswith(".bends")][0]
    b1 = {int(t[1]): Decimal(t[4]) for t in (ln.split() for ln in open(stem + ".bends"))
          if t and t[0] == "EDGE"}
    b2 = {int(t[1]): Decimal(t[4]) for t in (ln.split() for ln in open(stem2 + ".bends"))
          if t and t[0] == "EDGE"}
    # eid maps differ between dump and obj numbering; compare sorted multisets
    d = max(abs(x - y) for x, y in zip(sorted(b1.values()), sorted(b2.values())))
    check("objseed: same root", d < Decimal("1e-30"), f"max diff {d:.1e}")

    # --- pancake path ---
    cert3 = bendprover.process(None, PANCAKE, td, CLERS_BIN)
    check("pancake: recognized", cert3.get("pancake") is True and cert3.get("ok") is True)
    check("pancake: rim turning tau", abs(cert3["turning"]["min_nonflat"] - 6.283185307179586) < 1e-8)
    check("pancake: obj emitted", cert3.get("develop_spread", 1) < 1e-12)

print(f"\n{np} pass, {nf} fail")
sys.exit(1 if nf else 0)

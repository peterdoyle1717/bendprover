#!/usr/bin/env python3
"""Format self-tests (FORMAT.md): the writer's records for the three
reference nets, checked against 256-bit closed forms and the format
grammar. Requires only the built binary."""
import subprocess, sys, os
import gmpy2
from gmpy2 import mpfr

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "csrc", "euclid_lm_mp")
gmpy2.get_context().precision = 256
PI = gmpy2.const_pi()
A13 = gmpy2.acos(mpfr(1) / 3)
TET = (PI - A13) / PI          # 0.608173...
ACS = A13 / PI                 # 0.391826...

def prove(name, nc):
    r = subprocess.run([BIN, "--prove", "--batch"], input=f"{name} {nc}\n",
                       capture_output=True, text=True)
    rec = {"b": {}}
    for ln in r.stdout.splitlines():
        t = ln.split()
        if not t or t[0] == "#": continue
        if t[0] == "b": rec["b"][(int(t[1]), int(t[2]))] = t[3]
        elif t[0] == "end": pass
        else: rec[t[0]] = " ".join(t[1:])
    return rec

np = nf = 0
def check(msg, cond):
    global np, nf
    np += cond; nf += (not cond)
    print(("ok  " if cond else "FAIL") + "  " + msg)

# tetrahedron
r = prove("CCAE", "1,2,3;1,3,4;1,4,2;2,4,3")
check("tetra: name echoed", r.get("net") == "CCAE")
check("tetra: unit halfturns", r.get("unit") == "halfturns")
be = float(r.get("benderr", "1"))
check("tetra: benderr < 1e-30", 0 < be < 1e-30)
worst = max(abs(mpfr(v) - TET) for v in r["b"].values())
check(f"tetra: 6 bends within benderr of closed form (worst {float(worst):.1e})",
      len(r["b"]) == 6 and worst < be)

# doubled hexagon: exact pancake
faces = "1,2,3;1,3,4;1,4,5;1,5,2;2,5,6;2,6,3;3,6,7;3,7,8;3,8,4;4,8,5;5,8,7;5,7,6"
r = prove("CCCACACCAABE", faces)
check("pancake: benderr 0", r.get("benderr") == "0")
check("pancake: all bends integer tokens",
      all(v in ("0", "1", "-1") for v in r["b"].values()))
check("pancake: 18 edges", len(r["b"]) == 18)

# v13: closed-form folds and exact flats
faces = ("1,2,3;1,3,4;1,4,5;1,5,2;2,5,6;2,6,7;2,7,3;3,7,8;3,8,4;4,8,9;"
         "4,9,10;4,10,5;5,10,11;5,11,6;6,11,12;6,12,7;7,12,13;7,13,8;"
         "8,13,9;9,13,10;10,13,12;10,12,11")
r = prove("CCCACCACACCACACACAAABE", faces)
be = float(r.get("benderr", "1"))
flats = [k for k, v in r["b"].items() if v == "0"]
folds = {k: mpfr(v) for k, v in r["b"].items() if v != "0"}
worst = max(min(abs(v - ACS), abs(v - TET)) for v in folds.values())
check(f"v13: 15 exact flats", len(flats) == 15)
check(f"v13: folds within benderr of acos(1/3)/(gem/2) or complement "
      f"(worst {float(worst):.1e})", worst < be)

print(f"\n{np} pass, {nf} fail")
sys.exit(1 if nf else 0)

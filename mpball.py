"""mpball.py -- ball arithmetic for the embeddedness checker.

A Ball is (c, r): an mpfr center at BPREC bits and a float radius, with
the invariant that the represented exact value lies in [c-r, c+r]. Every
operation inflates r by a rigorous bound on the new uncertainty plus an
ulp term for the center rounding (2^(20-BPREC) * |c|, generous).

Only what the checker needs: +, -, *, neg, sqrt, sin, cos, division by
a positive ball, and 3-vector helpers.
"""
import math
import gmpy2
from gmpy2 import mpfr

BPREC = 192
gmpy2.get_context().precision = BPREC
ULP = 2.0 ** (20 - BPREC)

def _f(x):
    return float(x)

class Ball:
    __slots__ = ("c", "r")
    def __init__(self, c, r=0.0):
        self.c = mpfr(c); self.r = float(r)
    def __add__(a, b):
        b = b if isinstance(b, Ball) else Ball(b)
        c = a.c + b.c
        return Ball(c, a.r + b.r + ULP * (abs(_f(c)) + 1e-300))
    def __sub__(a, b):
        b = b if isinstance(b, Ball) else Ball(b)
        c = a.c - b.c
        return Ball(c, a.r + b.r + ULP * (abs(_f(c)) + 1e-300))
    def __neg__(a):
        return Ball(-a.c, a.r)
    def __mul__(a, b):
        b = b if isinstance(b, Ball) else Ball(b)
        c = a.c * b.c
        r = (abs(_f(a.c)) * b.r + abs(_f(b.c)) * a.r + a.r * b.r) * 1.0000001
        return Ball(c, r + ULP * (abs(_f(c)) + 1e-300))
    def divpos(a, b):
        """a / b where b is certified positive (b.c - b.r > 0)."""
        lo = _f(b.c) - b.r
        assert lo > 0, "divpos: divisor not positive"
        c = a.c / b.c
        r = (a.r + abs(_f(c)) * b.r) / lo * 1.0000001
        return Ball(c, r + ULP * (abs(_f(c)) + 1e-300))
    def sqrtpos(a):
        lo = _f(a.c) - a.r
        assert lo > 0, "sqrtpos: not positive"
        c = gmpy2.sqrt(a.c)
        r = a.r / (2.0 * math.sqrt(lo)) * 1.0000001
        return Ball(c, r + ULP * (abs(_f(c)) + 1e-300))
    def sin(a):
        return Ball(gmpy2.sin(a.c), a.r + ULP)
    def cos(a):
        return Ball(gmpy2.cos(a.c), a.r + ULP)
    def lo(a): return _f(a.c) - a.r
    def hi(a): return _f(a.c) + a.r
    def __repr__(a): return f"Ball({_f(a.c):.6g},{a.r:.2g})"

def bvec(x, y, z): return [Ball(x), Ball(y), Ball(z)]
def vadd(u, v): return [u[i] + v[i] for i in range(3)]
def vsub(u, v): return [u[i] - v[i] for i in range(3)]
def vscale(s, u): return [s * u[i] for i in range(3)]
def vdot(u, v): return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def vcross(u, v):
    return [u[1]*v[2] - u[2]*v[1],
            u[2]*v[0] - u[0]*v[2],
            u[0]*v[1] - u[1]*v[0]]
def vnormpos(u):
    return vdot(u, u).sqrtpos()
def vunit(u):
    n = vnormpos(u)
    return [u[i].divpos(n) for i in range(3)]
def vmaxrad(u): return max(x.r for x in u)

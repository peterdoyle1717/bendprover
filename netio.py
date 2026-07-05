"""netio.py -- net loading, seed bends, development, output writers.

Conventions:
- Nets are unit-equilateral triangulated spheres. Faces are oriented ccw
  viewed from outside. Edge ids are canonical sorted vertex pairs.
- Bend = pi - dihedral, signed; 0 = flat, tau/2 = folded shut (pancakes only).
- Vertex stars are cyclic edge lists consistent with the face orientation.
- Closure model (validated at 1e-15 on the 322-net census):
      H(b) = prod over the star of Rz(60deg) Rx(b_e) = Id.
"""
import os, math, subprocess

# ---------- CLERS ----------
def faces_from_clers(clers, clers_bin="clers"):
    """decode via the clers binary (reads the name on stdin,
    writes '1,2,3;1,3,4;...' -- 1-based, orientation-consistent)."""
    out = subprocess.run([clers_bin, "decode"], input=clers + "\n",
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"clers decode failed: {out.stderr.strip()[:120]}")
    faces = []
    for chunk in out.stdout.strip().replace("\n", ";").split(";"):
        chunk = chunk.strip()
        if not chunk: continue
        parts = [int(x) for x in chunk.replace(",", " ").split()]
        for i in range(0, len(parts), 3):
            faces.append((parts[i], parts[i+1], parts[i+2]))
    return faces

# ---------- combinatorics ----------
def stars_from_faces(faces):
    """eab: eid -> (a,b) sorted; stars: vid -> cyclic edge list (orientation-consistent)."""
    apex = {}
    for a, b, c in faces:
        apex[(a, b)] = c; apex[(b, c)] = a; apex[(c, a)] = b
    und = sorted(set(tuple(sorted(k)) for k in apex))
    eid = {p: i for i, p in enumerate(und)}
    eab = {i: p for p, i in eid.items()}
    verts = sorted(set(v for f in faces for v in f))
    stars = {}
    for v in verts:
        nbrs = [w for (x, w) in apex if x == v]
        w0 = nbrs[0]
        cyc = [w0]
        while True:
            w = apex[(v, cyc[-1])]
            if w == w0: break
            cyc.append(w)
        stars[v] = [eid[tuple(sorted((v, w)))] for w in cyc]
    return eab, stars

# ---------- dump seeds (bend_sweep format) ----------
def read_dump(path):
    """returns clers, eab, bends(str), stars -- stars from the dump's VERT lines."""
    clers = None; eab = {}; bends = {}; stars = {}
    for ln in open(path):
        t = ln.split()
        if not t: continue
        if t[0] == "CASE":
            clers = t[1]
        elif t[0] == "EDGE":
            eab[int(t[1])] = (int(t[2]), int(t[3]))
            bends[int(t[1])] = t[4]
        elif t[0] == "VERT":
            stars[int(t[1])] = [int(x) for x in " ".join(t[2:]).split("|")[0].split(",") if x.strip()]
    return clers, eab, bends, stars

# ---------- obj seeds ----------
def read_obj(path):
    V = []; F = []
    for ln in open(path):
        t = ln.split()
        if not t: continue
        if t[0] == "v": V.append(tuple(float(x) for x in t[1:4]))
        elif t[0] == "f": F.append(tuple(int(x.split("/")[0]) - 1 for x in t[1:4]))
    return V, F

def _sub(u, v): return (u[0]-v[0], u[1]-v[1], u[2]-v[2])
def cross(u, v): return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
def dot(u, v): return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]
def norm(u): return math.sqrt(dot(u, u))

def bends_from_obj(path):
    """signed bends from coordinates; sign fixed globally by closure test upstream."""
    V, F = read_obj(path)
    faces = [tuple(f) for f in F]
    eab, stars = stars_from_faces([(a+1, b+1, c+1) for a, b, c in faces])  # 1-based vertices
    # face normals
    nrm = {}
    for f in faces:
        a, b, c = (V[f[0]], V[f[1]], V[f[2]])
        n = cross(_sub(b, a), _sub(c, a))
        nn = norm(n)
        nrm[f] = (n[0]/nn, n[1]/nn, n[2]/nn)
    side = {}
    for f in faces:
        a, b, c = f
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            side[(x+1, y+1)] = f
    bends = {}
    for e, (p, q) in eab.items():
        f1 = side.get((p, q)); f2 = side.get((q, p))
        n1, n2 = nrm[f1], nrm[f2]
        cosd = max(-1.0, min(1.0, dot(n1, n2)))
        edge = _sub(V[q-1], V[p-1])
        s = dot(cross(n1, n2), edge)
        b = math.atan2(math.copysign(norm(cross(n1, n2)), s), cosd)
        bends[e] = repr(b)
    return eab, bends, stars

# ---------- development (bends -> coords) ----------
S3 = math.sqrt(3.0)/2
def _scl(s, u): return (s*u[0], s*u[1], s*u[2])
def _add(u, v): return (u[0]+v[0], u[1]+v[1], u[2]+v[2])
def _unit(u): n = norm(u); return (u[0]/n, u[1]/n, u[2]/n)

def _rot_about(pt, ax, th, p):
    u = _unit(ax); v = _sub(p, pt)
    c, s = math.cos(th), math.sin(th)
    r = _add(_add(_scl(c, v), _scl(s, cross(u, v))), _scl(dot(u, v)*(1-c), u))
    return _add(pt, r)

def develop(faces, eab, bends_f, sign=+1.0):
    """coords per vertex from bends; returns (coords, worst vertex spread).

    Unfolds face-by-face across DIRECTED edges (the placed face's winding),
    so the rotation axis direction is per-edge consistent and a single
    global sign covers the convex/reflex convention (develop_best tries both).
    """
    from collections import deque
    apex = {}
    for a, b, c in faces:
        apex[(a, b)] = c; apex[(b, c)] = a; apex[(c, a)] = b
    pair2e = {tuple(sorted(v)): k for k, v in eab.items()}
    fkey = lambda f: min(f, f[1:] + f[:1], f[2:] + f[:2])
    a, b, c = faces[0]
    first = {a: (0.0, 0.0, 0.0), b: (1.0, 0.0, 0.0), c: (0.5, S3, 0.0)}
    fpos = {fkey(faces[0]): (faces[0], first)}
    dq = deque([faces[0]])
    while dq:
        f0 = dq.popleft()
        mp = fpos[fkey(f0)][1]
        for i in range(3):
            x, y = f0[i], f0[(i + 1) % 3]
            d = apex[(y, x)]
            f1 = (y, x, d)
            if fkey(f1) in fpos: continue
            A, Bp, C1 = mp[x], mp[y], mp[f0[(i + 2) % 3]]
            u = _unit(_sub(Bp, A))
            M = _add(A, _scl(dot(_sub(C1, A), u), u))
            D0 = _sub(_scl(2, M), C1)
            bend = bends_f[pair2e[tuple(sorted((x, y)))]]
            D = _rot_about(A, u, sign * bend, D0)
            fpos[fkey(f1)] = (f1, {y: Bp, x: A, d: D})
            dq.append(f1)
    from collections import defaultdict as dd
    pos = dd(list)
    for f, mp in fpos.values():
        for vid, p in mp.items(): pos[vid].append(p)
    worst = 0.0; coords = {}
    for vid, ps in pos.items():
        cx = tuple(sum(p[i] for p in ps) / len(ps) for i in range(3))
        coords[vid] = cx
        for p in ps: worst = max(worst, norm(_sub(p, cx)))
    return coords, worst

def develop_best(faces, eab, bends_f):
    best = None
    for sign in (+1.0, -1.0):
        coords, worst = develop(faces, eab, bends_f, sign)
        if best is None or worst < best[1]: best = (coords, worst, sign)
    return best

# ---------- writers ----------
def write_obj(path, faces, coords):
    verts = sorted(coords)
    idx = {v: i+1 for i, v in enumerate(verts)}
    with open(path, "w") as f:
        for v in verts:
            f.write("v %.17g %.17g %.17g\n" % coords[v])
        for a, b, c in faces:
            f.write(f"f {idx[a]} {idx[b]} {idx[c]}\n")

def write_bends(path, clers, eab, flats, fold_bends, cert, turning):
    """fold_bends: eid -> 40-digit string; flats written as exact 0."""
    with open(path, "w") as f:
        f.write(f"CLERS {clers}\n")
        f.write(f"CERT {'PASS' if cert.get('ok') else cert.get('why','FAIL')} "
                f"sigma_low {cert.get('sig_low','-')} h {cert.get('h','-')} radius {cert.get('radius','-')}\n")
        f.write(f"TURNING min_nonflat {turning.get('min_nonflat','-')} flats_exact_zero {turning.get('nflatv',0)}\n")
        for e in sorted(eab):
            a, b = eab[e]
            if e in flats:
                f.write(f"EDGE {e} {a} {b} 0\n")
            else:
                f.write(f"EDGE {e} {a} {b} {fold_bends[e]}\n")

# bendprover

Solver + certificate for the bends of unit-equilateral triangulated
spheres ("nets"), at uniform 128-bit precision.

Given a net and a numerical seed, it classifies edges as flat or
folded, solves the vertex-closure system for the fold bends, and
certifies existence and local uniqueness of a nearby exact root.
Flat edges are reported as exactly 0.

## Formulation

Fold edge `e` carries a quaternion pair `(c,s) ~ (cos b/2, sin b/2)`
with unit row `c² + s² = 1`. Vertex closure is the quaternion word

    ∏ over the star of  T₆₀ · (c,s,0,0)  =  ±1,

with `T₆₀ = (√3/2, 0, 0, 1/2)` the 60° face turn and flat edges
structurally absent (faces merged into flat plates). The residuals are
polynomial — multilinear in the pairs — so the prover needs no
transcendental functions. Newton runs at 128 bits (gmpy2/MPFR).

The certificate takes a numerically selected square row subsystem,
lower-bounds its smallest singular value by a Cholesky factor witness
(`δ = ‖I−BC‖`, `σ ≥ (1−δ)/‖B‖`, margin `σ_C² − ‖A−CCᵀ‖`), and applies
Kantorovich with the explicit multilinear Hessian bound `L = 4·arity²`.
Rounding is charged to the standard 128-bit model with a ×64 slop
factor. What is certified: an exact root of the selected subsystem
within the reported radius; the dropped closure rows are evaluated and
bounded (`drop_bound`), and their exact discharge is the lune lemma
plus the monodromy identity at one excised plate.

Pancakes — nets whose solution is a doubly covered flat region, every
edge 0 or τ/2 — are recognized and passed through without solving;
they are included in the outputs.

## Requirements

python3, gmpy2 (`pip install gmpy2`), and the `clers` binary
(CLERS name → face list) on PATH or named by `CLERS_BIN`.

## Usage

    # single net from a bend_sweep dump (CLERS read from its CASE line)
    python3 bendprover.py --dump net.dump OUT/

    # single net seeded from coordinates
    python3 bendprover.py --clers CCC...E --obj net.obj OUT/

    # batch: one seed path per line (or "CLERS<TAB>path")
    python3 bendprover.py --batch list.txt OUT/

Per net `OUT/` gets `v<NV><CLERS>.bends` (flats exact 0, folds to
40 digits, certificate and turning summary in the header),
`.cert` (full certificate record, JSON), and `.obj` (coordinates
developed from the bends, for the embeddedness step; per-vertex
placement spread is reported as `develop_spread`).

The `TURNING` line reports the minimum total turning over non-flat
vertices; positive means the object is neoconvex with exact flats.

Obj seeds must use the same vertex numbering as the `clers` decode of
the named net (the pipeline checks and refuses otherwise); dumps carry
their own numbering and always work.

## Tests

    make test

runs an end-to-end fold net (certify, bends match a 200-digit
reference dump, obj round-trip re-certifies the same root, turning
positive) and the pancake path.

## Reference run

`census322/` holds the output of a full run over the 322-net flopper
census (276 certified fold nets + 46 pancakes), produced by
`make census` with the dump directory named in the Makefile.

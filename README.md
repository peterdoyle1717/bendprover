# bendprover

Solver and existence prover for the bends of unit-equilateral
triangulated spheres ("nets"), at uniform multiprecision.

The pipeline is one C program, [csrc/euclid_lm_mp.c](csrc/euclid_lm_mp.c):
Levenberg–Marquardt with a dent gate on the quaternion holonomy system,
from a closed-form wish start, at `--prec` bits (default 128). In
`--prove` mode each net gets the full treatment in one invocation:

    solve from scratch -> classify flat edges (1e-8) -> re-solve with
    flats frozen to exactly 0 -> take-1 existence certificate
    (Newton polish on quaternion pairs, numerically selected square
    subsystem, Cholesky factor witness, Kantorovich) -> emit record

Records follow [FORMAT.md](FORMAT.md): bends in halfturns, integer
values exact (0 = flat, 1 = folded shut; degenerate pancakes are
all-integer with `benderr 0`), decimals certified within the record's
`benderr`, proof internals in a `#` provenance comment.

## Build and run

    make euclid_lm_mp          # needs mpfr + gmp
    echo 'CCAE 1,2,3;1,3,4;1,4,2;2,4,3' | csrc/euclid_lm_mp --prove --batch

Batch lines are `NAME NETCODE` (name optional). `run/doob_prove.sh`
is the production driver: chunked netcode lists, one solver process
per chunk, nice 19.

## Downstream checks

Records are the sole deliverable; consumers need nothing else.

- neoconvexity: the turning window per vertex is reported in the
  provenance comment and rechecked trivially from the bends.
- embeddedness: [csrc/embcheck_mp.c](csrc/embcheck_mp.c) develops
  ball-arithmetic vertex enclosures from (bends, benderr) and
  certifies shared-edge, vertex-link, and pairwise-separation
  conditions. Verdicts per net: `PASS` (embedded certified),
  `PANCAKE` (all bends integer halfturns; the flat doubly-covered
  case, certified by exactness, not embedded), `FAIL` (nothing
  certified). The supporting lemma (reduction of shared-vertex
  crossings) is documented in `embcheck.py`'s docstring and still
  needs a written proof. `run/doob_emb.sh` sweeps a prove run's
  chunks.
- objs: coordinates are derived data; `euclid_lm_mp NETCODE` without
  `--bends-only` develops and prints them on demand.

## Python oracle

`mpkernel.py`, `netio.py`, `bendprover.py`, `embcheck.py`/`mpball.py`,
and `tests/` are the reference implementations the C programs were
ported from. They run nothing in production; `make test-oracle`
exercises them, and the C binaries are differential-tested against
them (embcheck: 0 verdict mismatches on 365 mixed census/sweep nets).

## Tests

    make test          # format self-tests: closed-form bend checks
    make test-oracle   # the Python reference pipeline

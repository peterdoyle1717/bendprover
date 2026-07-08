# Bend record format

One record per net. Line-oriented, self-delimiting, designed to be
readable without this document.

```
net CCAE
v 4
e 6
unit halfturns
benderr 1.8e-34
faces 1,2,3;1,3,4;1,4,2;2,4,3
# proof: kantorovich at 128 bits, jacobian sigma_min 3.09e-01, h 2.3e-31, solver euclid_lm_mp
b 1 2 0.608173447969392729829144440777569089
b 1 3 0.608173447969392729829144440777569089
b 1 4 0.608173447969392729829144440777569089
b 2 3 0.608173447969392729829144440777569089
b 2 4 0.608173447969392729829144440777569089
b 3 4 0.608173447969392729829144440777569089
end
```

## Grammar

A record is:

- `net NAME` — canonical CLERS name (program-wide identity);
- key-value lines `v N`, `e N`, `unit halfturns`, `benderr X`;
- `faces a,b,c;d,e,f;...` — the labeled triangulation. THIS is the
  labeling the `b` lines refer to; no decoder is needed downstream.
  Faces are wound counterclockwise viewed from outside;
- `#` lines — provenance comments, no semantics;
- `b A B VALUE` lines, one per edge, sorted by vertex pair (A < B,
  lexicographic);
- `end`.

Records concatenate freely into stream files. Failed nets appear as a
single `# failed NETCODE reason` comment line.

## Semantics

- Bends are in halfturns: VALUE = (signed exterior dihedral)/(gem/2).
  Positive = convex, under the ccw-outside orientation. Range (-1, 1].
- **Integer values are exact claims**: `0` = the edge is exactly flat,
  `1` (or `-1`) = folded exactly shut. They are not approximations.
- **Decimal values are certified**: each is within `benderr` (absolute,
  in halfturns) of the corresponding bend of an exact unit-equilateral
  realization of the net. `benderr 0` means every value in the record
  is exact; such records contain only integer values (the degenerate
  pancakes).
- benderr is derived from the existence certificate's ball radius r in
  pair space: |bend error| <= 2r/(1-r) in radians (the gradient of
  2*atan2 on the certified ball, where the unit rows keep
  c^2+s^2 >= (1-r)^2), divided by gem/2, plus the print quantum
  0.5e-36. The proof internals live in the `#` provenance comment.

## Self-tests

Any writer must reproduce, within its stated benderr:

- tetrahedron (CCAE): all six bends (gem/2 - acos(1/3))/(gem/2)
  = 0.608173447969392729829144440777569089...
- the doubled hexagon (CCCACACCAABE): benderr 0, rim bends exactly 1,
  interior exactly 0;
- the v13 net CCCACCACACCACACACAAABE: fold values acos(1/3)/(gem/2)
  = 0.391826552030607270170855559222430911... and its complement
  0.608173..., flats exactly 0.

## Manifest

Stream files are registered NEO_DATA-style, one line each:

```
file=NAME.bends  records=N  max_benderr=X  sha256=H  date=YYYY-MM-DD
```

The manifest is an audit aggregate; every bound a consumer needs is in
the records themselves.

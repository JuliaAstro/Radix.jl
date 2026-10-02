# Field names for Radix rate types

These names are applied in `src/`. This file records the rules, the shared
vocabulary, what was checked against the XSTAR Fortran, and the resulting
fields of every type (the per-type tables are generated from the structs).

The original names came from the XSTAR data-type tables (`i1…i8`, `r1…`,
`ionN`, `kN−1`). Those are positional labels and were inconsistent across types:
`N` meant ion index in one type and electron count in another, and the same
level was `i`, `k`, `j`, `n` or `iN`.

## Rules

1. **Name the meaning, not the Fortran slot.** Words for roles; single letters
   only for standard physics symbols (`Z`, `n`, `A`, `λ`, `σ`, `E`, `T`).
2. **One name per concept** across all types (vocabulary below).
3. **Every record starts with `rtype::Int8` and `label::String`.** `rtype` is
   XSTAR's rate type (`lrtyp`, e.g. 5, 8, 13). It is not the data-type code and
   it matters (`ChargeExH0` branches on it). It used to be called `rate`, which
   shadowed the function `rate`.
4. **Concrete, typed containers.** Tabulated data and fit coefficients are
   `Vector{R}`, two-dimensional tables are `Matrix{R}`; no `Tuple`. Sizes
   (`nd`, `nt`, `nx`) are not stored, they are the array lengths. Integers are
   `I<:Integer`, floats `R<:AbstractFloat`.
5. **Units in comments, not names.** Energies eV unless noted (`ΔE` in Ry for
   the CHIANTI fits). Temperatures K, cross sections Mb, rates s⁻¹ or cm³ s⁻¹.
   `Cell.T` is in 10⁴ K as in XSTAR.

   `Cell` holds the local gas state only: `T`, `nₕ` (`xh0`), `nₑ` (`xnx`) and
   `ntot` (`xpx`, the hydrogen density). Level data and element masses live in
   tables (`level_table`), and the radiation field will be a separate object.
6. **Greek only for standard, easy-to-type symbols** (`λ σ α γ ρ η Υ`). Only
   `E∞` is spelled out (`E_inf`), because a subscript is not a valid identifier
   character.
7. **Group values that travel together in small structs** (below).

## Shared vocabulary

| Concept | Name | Replaces |
|---|---|---|
| atomic number | `Z` | `Z` |
| ion index (global, see below) | `ion` | `N`, `ionN`, `N_`, `NN` |
| single level index | `level` | `i`, `n` (AtomicLevel), `iN` |
| transition levels | `transition::Transition(lower, upper)` | `i`/`k`, `i`/`j` |
| parent ion and level | `parent::Parent(ion, level)` | `Nm`, `ionNm`, `kNm`, `km`, `k` |
| electrons of the recombining ion (N−1) | `parent_electrons` | `Nm` (TotRadRecomb) |
| electrons in the ion | `n_electrons` | `NN` |
| principal quantum number | `n` | `n` |
| orbital angular momentum (term) | `L` | `L` |
| orbital quantum number of a subshell | `l` | `L` (ParPhotoIonize3) |
| spin multiplicity 2S+1 | `spin_mult` | `S2`, `S2p` |
| 2J | `twoJ` | `J2` |
| statistical weight 2J+1 | `g` | `J2p`, `J2` (level data) |
| effective quantum number | `n_eff` | `ν` (misdocumented as frequency) |
| transition type flag | `kind` | `type`, `t`, `it` |
| superlevel index (0 if absent) | `superlevel` | `k` (PhotoionizeDamp) |
| energy | `E` | `E`, `Ei`, `e` |
| threshold / ionization / limit energy | `E_th`, `E_ion`, `E_inf` | `Eth`, `Ei`, `Einf`, `E∞` |
| transition energy | `ΔE` | `ΔE` |
| wavelength | `λ` | `λ`, `e` |
| Einstein A | `A` | `A`, `a` |
| autoionization / radiative rates | `A_auto`, `A_rad` | `Aa`, `Ar` |
| (weighted) oscillator strength | `gf` / `f` | `gf`, `f`, `g` |
| tabulated grids | `T_grid` (log₁₀ K for `ElectronImpact1`), `ne_grid`, `E_grid` | `T`, `Te`, `ne`, `E` (tuples) |
| cross section | `σ` | `σ` |
| fit coefficients (unnamed) | `coeffs` | `c::Tuple` |
| effective collision strength | `Υ` (U+03A5) | `Υ` |
| named fit parameters | `a b c d e`, `A B C`, `T0 T1 T2` | same |

## Transition and Parent

`Transition(lower, upper)` and `Parent(ion, level)` live in
`src/transition.jl`. A line or collision record is `coef.transition.lower` /
`.upper`; a photoionization or Auger record is `coef.parent.ion` / `.level`.

`Parent.ion` is `0` when the data do not give the parent ion. This applies to
`ChargeExHe`, `PhotoionizeDamp` and `AutoionizeSat`; they use `Parent` anyway so
the code is the same everywhere and only `parent.level` is read.
`TotRadRecomb` is the one exception: it stores `parent_electrons`, an electron
count, not an ion index.

## Verified against the Fortran and the data

Checked against `ucalc()` in `xstarsub.f` and by loading `atdb.fits`:

- **`ion` is a global ion index** in the ion table (h_i=1, he_i=2, he_ii=3,
  li_i=4 … fe_xxvi=351), not an electron count. `Ion.N` was this index and
  `Ion.I` is the stage (1 = neutral), so they are now `ion` and `stage`.
- **The parent ion is `ion + 1`** for 287,173 of 287,180 `ParPhotoIonize1`
  records (the remainder are special cases such as fully stripped ions), so it
  is kept as data, not derived.
- **`TotRadRecomb.Nm` is an electron count (N−1)**, not an ion index →
  `parent_electrons`. `ParPhotoIonize3.NN` is the electron count of the ion →
  `n_electrons`.
- **`J2p` in `AtomicLevel`/`AtomicLevelFe` is 2J+1** (2P1/2 → 2, 2P3/2 → 4), the
  statistical weight → `g`. Other types store 2J as an integer → `twoJ`.
- **`rate::Int8` is the record's rate type (`lrtyp`)**, not the data-type code
  → `rtype`; kept because code uses it.
- **Charge exchange.** Type 9 (`ChargeExHe`) takes `idat(1)` as its own level
  and `idat(2)` → `nlevp + idat(2) − 1`, a level of the recombined ion; with one
  integer both default to 1. Type 10 (`ChargeExHp`) ignores its integers in
  `ucalc`; in the data they are (level, ion index of the neutral, e.g. 29 =
  O I). Both are inferred.
- **`PhotoionizeDamp.k`** is a superlevel index present only in the 8-integer
  variant and 0 otherwise → `superlevel`.
- **`Atom`:** both integers equal `Z` in every record (`n_ions`, `Z`); the first
  float is the abundance (H 1.0, He 0.1, C 3.7e-4), the second the atomic mass.
- **`TotRadRecombH`:** the Fortran uses the first integer, `Z`, as the effective
  charge, so `rate` uses `coef.Z^2`.
- **Reading `atdb.fits`.** FITSFiles.jl up to v0.3.2 rounds unscaled Int32
  columns through Float32, so the pointers to the reals were wrong beyond 2^24
  for about 1M of the 1.2M records (reals read from the wrong place, e.g. S II
  levels with `g = 0`). `load` passes `scale=false` (the file has no
  TSCAL/TZERO); fixed upstream in JuliaAstro/FITSFiles.jl#50. Values read from
  records past #71,119 before this change were wrong.
- **Equivalence of the refactor.** Loading `atdb.fits` before and after
  introducing `Transition`/`Parent` (and moving the three no-parent-ion types
  onto `Parent`) gave identical lower, upper, parent-ion and parent-level values
  for every record of all 45 types that hold data.

## Bugs found and fixed while renaming

- `ratemap[9]` pointed at `ChargeExH0` instead of `ChargeExHe` (59 records).
- `CollisionHelikeSat` records were built as `CollisionHelike` (427 records).
- `TwoPhotonRad` called an undefined `TwoPhoton`.
- `TotRadRecombH` used an undefined `nmax` (now `Z`).

## Open items

- `CollisionIonize` stores `ρ` from `rvec[2:end]`; the header says it starts at
  `r3`.
- `ChargeExH0.rate` clamps with `max(0, …)` and tests `T > 5`; neither is in
  `ucalc`.
- The `Atom` records for iron and zinc look wrong (abundance 3.48 and 9.32, mass
  3.70 and 1.0); check the source table.
- `AtomicLine.rate` is unfinished and only had its field names updated.
- The integer meanings of `ChargeExHe` and `ChargeExHp` are inferred.
- `DielecRecombH` and `TotDielecRecomb` store no data yet.

## Fields of every type

`rtype` and `label` come first in every type and are omitted. Numbers in
brackets are the XSTAR data-type codes.

**Levels and lines**

| Type (code) | Fields |
|---|---|
| `AtomicLevel` (6) | `n`, `spin_mult`, `L`, `Z`, `level`, `ion`, `E`, `g`, `n_eff`, `E_inf` |
| `AtomicLevelFe` (83) | `level`, `ion`, `E`, `g` |
| `AtomicLine` (4) | `transition`, `λ`, `f`, `A` |
| `AtomicLine2` (50) | `transition`, `Z`, `ion`, `λ`, `gf`, `A` |
| `RadiativeAPED` (91) | `transition`, `Z`, `ion`, `λ`, `A` |
| `RadiativeProb` (54) | `transition`, `Z`, `ion`, `A` |
| `RadiativeFeDecay` (82) | `transition`, `ion`, `λ`, `E`, `gf`, `A_rad`, `A_auto` |
| `TwoPhotonDecay` (76) | `transition`, `ion`, `A` |
| `TwoPhotonRad` (5) | `transition`, `r1`, `r2` |
| `EffectiveCharge` (57) | `n`, `L`, `twoJ`, `Z`, `level`, `ion`, `Zeff` |

**Autoionization**

| Type (code) | Fields |
|---|---|
| `Autoionize` (3) | `C`, `E` |
| `AutoionizeSat` (72) | `spin_mult`, `L`, `level`, `parent`, `Z`, `ion`, `A_auto`, `E`, `g` |
| `AutoionizeFe25Sat` (75) | `ion`, `level`, `parent`, `A_auto`, `E` |
| `IronKAuger` (86) | `parent`, `level`, `Z`, `ion`, `E`, `A_widths` |

**Recombination and charge exchange**

| Type (code) | Fields |
|---|---|
| `RadRecomb` (1) | `A`, `η` |
| `TotRadRecomb` (38) | `Z`, `parent_electrons`, `M`, `W`, `ion`, `A`, `B`, `T0`, `T1`, `C`, `T2` |
| `TotRadRecombH` (30) | `Z`, `ion` |
| `DielecRecomb1` (7) | `A`, `B`, `T0`, `T1` |
| `DielecRecomb2` (8) | `C`, `E` |
| `DielecRecombH` (22) | `rate` |
| `TotDielecRecomb` (39) | `rate` |
| `ChargeExH0` (2) | `ion`, `a`, `b`, `c`, `d`, `T1`, `T2`, `ΔE` |
| `ChargeExHe` (9) | `level`, `parent`, `a`, `b`, `c`, `d`, `T1`, `T2`, `ΔE` |
| `ChargeExHp` (10) | `level`, `ion`, `a`, `b`, `c`, `d`, `e`, `T1`, `T2`, `ΔE` |

**Collisions**

| Type (code) | Fields |
|---|---|
| `CollisionHlike1` (60) | `transition`, `ion`, `coeffs` |
| `CollisionHlike2` (62) | `transition`, `ion`, `coeffs` |
| `CollisionHelike` (68) | `transition`, `Z`, `ion`, `coeffs` |
| `CollisionHeFine` (66) | `transition`, `Z`, `ion`, `coeffs` |
| `CollisionHelikeSat` (73) | `transition`, `Z`, `ion`, `coeffs` |
| `CollisionLS` (69) | `transition`, `Z`, `ion`, `coeffs` |
| `CollisionFe19` (81) | `transition`, `Z`, `ion`, `Υ` |
| `CollisionProb` (63) | `transition`, `Z`, `ion` |
| `CollisionAPED` (92) | `transition`, `Z`, `ion`, `T_grid`, `Υ` |
| `CollisionIonize` (95) | `level`, `ion`, `E_th`, `T0`, `ρ` |
| `ElectronCollision` (51) | `kind`, `transition`, `Z`, `ion`, `ΔE`, `C`, `Υ` |
| `ElectronImpact1` (56) | `transition`, `Z`, `ion`, `T_grid`, `Υ` |
| `ElectronImpact2` (98) | `kind`, `transition`, `ion`, `ΔE`, `C`, `Υ` |

**Photoionization**

| Type (code) | Fields |
|---|---|
| `ParPhotoIonize1` (49) | `n`, `L`, `twoJ`, `Z`, `parent`, `level`, `ion`, `E_grid`, `σ` |
| `ParPhotoIonize2` (53) | `n`, `L`, `twoJ`, `Z`, `parent`, `level`, `ion`, `E_grid`, `σ` |
| `ParPhotoIonize3` (59) | `n_electrons`, `n`, `l`, `parent`, `level`, `ion`, `E_th`, `E0`, `σ0`, `ya`, `P`, `yw` |
| `PhotoionizeFeKedge` (85) | `n`, `L`, `twoJ`, `Z`, `parent`, `level`, `ion`, `Zeff`, `E_th`, `f`, `γ`, `scale` |
| `PhotoionizeDamp` (88) | `n`, `L`, `twoJ`, `Z`, `parent`, `superlevel`, `level`, `ion`, `E_grid`, `σ` |
| `PhotoionizeDelta` (74) | `n`, `L`, `spin_mult`, `Z`, `parent`, `level`, `ion`, `E_inf`, `E_grid`, `f` |

**Superlevels (tables over `ne_grid × T_grid`)**

| Type (code) | Fields |
|---|---|
| `RadiativeSuper` (71) | `transition`, `Z`, `ion`, `λ`, `ne_grid`, `T_grid`, `A` |
| `CollisionSuper` (77) | `transition`, `Z`, `ion`, `λ`, `ne_grid`, `T_grid`, `C` |
| `PhotoionizeSuper` (70) | `n`, `L`, `spin_mult`, `Z`, `parent`, `level`, `ion`, `ne_grid`, `T_grid`, `logα`, `E_grid`, `σ` |
| `PhotoRecombX` (99) | `n`, `L`, `spin_mult`, `Z`, `parent`, `level`, `ion`, `ne_grid`, `T_grid`, `α`, `E_grid`, `σ` |

**Atom and Ion**

| Type | Fields |
|---|---|
| `Atom` (13) | `n_ions, Z, abundance, mass` |
| `Ion` (14) | `stage, Z, ion, E_ion` |

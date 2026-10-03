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
8. **No numeric literals in expressions.** Every empirical or physical number is a
   named constant. Numbers used by more than one file are in `src/constants.jl`
   (`Ry_eV`, `hc_eVÅ`, `T_unit`, `kT_eV`, `collision_rate_coeff`, `cx_unit`, `tiny`, …);
   numbers that belong to one rate or fit are `const`s at the top of its file
   (`cxHe_fraction`, `ps_rho_coeff`, `sz_abethe`, …). Only structural numbers stay
   (0, 1, 2, powers, halving, array indices).

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

Checked against `ucalc.f90` (the Fortran 90 in `ftools/xstar/xstarlib/src`, which is what
the `xstar` executable is built from; `utils/xstarsub.f` is an older copy that is not
built and differs in constants and logic) and by loading `atdb.fits`:

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

## Bugs found and fixed

- `ratemap[9]` pointed at `ChargeExH0` instead of `ChargeExHe` (59 records).
- `CollisionHelikeSat` records were built as `CollisionHelike` (427 records).
- `TwoPhotonRad` called an undefined `TwoPhoton`.
- `TotRadRecombH` used an undefined `nmax` (now `Z`).
- Found by comparing with the real `ucalc`: `ChargeExHe` left out the factor `0.1 nₕ`
  (and put the rate in the wrong direction); `TotRadRecomb` used its characteristic
  temperatures in K against a temperature in 10⁴ K; `ElectronImpact2` read its
  integers as `[it, i, k, ionN]` when the record is `[i, k, it, ionN]` and took the
  scale `C` from the wrong real.
- The ports made from the older `xstarsub.f` used the wrong constants (13.598 eV per Ry,
  hc = 12398.54), a wrong sign in the third segment of the 5-point spline, and
  an out-of-date `AtomicLine2`; all aligned with `ucalc.f90`.

## Radiation and Opacity

Photoionization needs the radiation field, which is not part of `Cell`:
`Radiation(E, F)` holds XSTAR's energy grid (`epi`, eV, geometrically spaced) and the
spectrum on it (`bremsa`). The rates take it as `radiation=`, the number of levels of
the ion as `nlev=` (`level_counts(levels)` gives it for every ion), the escape
probabilities as `ptmp=` and the populations as `abund=`. `Opacity(n)` holds the
continuum opacity and recombination emissivity arrays that the integrals add to when
passed as `opacity=`. The shared engine is `photoionize_level` (types 49 and 53) with
`photoionization_integrals` (XSTAR's `phint53`) and `phextrap`.

## Reference tests

`test/ucalc_tests.jl` compares each ported rate with XSTAR's real `ucalc` on records
sampled from `atdb.fits` (`test/reference/ucalc/`, with the driver, the build script
and the generator). The driver is built like HEASoft's `xstar`, so XSTAR's
single-precision literals limit the agreement to about 1e-6. Types covered: 1, 2, 9, 30,
38, 49, 50, 51 (5- and 9-point fits), 53, 54, 56, 57, 59, 63, 70, 71, 74, 85, 86, 88, 98 and 99.

### Superlevels (types 70 and 99)

`PhotoionizeSuper` (70) and `PhotoRecombX` (99) tabulate a recombination coefficient on a
(log₁₀ nₑ, log₁₀ T) grid and a cross-section table. `rate` interpolates the coefficient,
scales the table so that the Milne relation (`milne_recombination`) reproduces it, integrates it
with `photoionization_integrals_hunt` (XSTAR's `phint53hunt`) and scales the photoionization rate
and heating by the ratio of the tabulated to the integrated recombination (`photoionize_superlevel`).
Both use the last real level (`nlev − 1`) at most as the lower level. Differences and quirks, all
reproduced:

- Type 70 reads its coefficient table with the temperature running fastest (`logα` has size
  `(nT, nₑ)`) and type 99 with the density running fastest (`α` has size `(nₑ, nT)`); the type 70
  code limits its second tabulated density to 10⁸ and, with `neutral=true` (`jkion = 1` in ucalc),
  the density to 10⁸ cm⁻³; it caps the scaled cross section at 10⁶ Mb and drops the tail below 10⁻⁶
  of its first point.
- Type 99 works in single precision. It interpolates in log₁₀ α linearly in T and, only for the
  density intervals 2 to nₑ−1, linearly in log₁₀ nₑ; at the lowest density, below it and above the
  table it uses the first density column. log₁₀ T is limited to the table.
- Neither scales the recombination cooling (`ienergy`) with the other rates, and an energy
  correction is applied to `fenergy2` and `ienergy2` in type 99 only.
- When the integrated recombination is below 10⁻⁴⁸ ucalc returns no rates but leaves the
  integrals `piht2` and `rrcl2` in its energy outputs `ans5` and `ans6`; `ienergy2` and `fenergy2`
  show them (with the sign of those outputs).
- `phint53hunt` reuses the bins of earlier passes, but not the weight `atmp22` of `rrcl2`, which keeps
  the value of the last bin computed in the pass, and with a span of one bin it integrates nothing.
- All the superlevel records of `atdb.fits` leave the ground level of the parent ion; the
  excited-parent branch is compared with a made-up parent level (`type99par`).

## Automatic differentiation

The rates are generic in the number type of the `Cell`, so ForwardDiff can differentiate
`frate` and `irate` with respect to the temperature or the densities: build the `Cell` from dual
numbers (all four fields of the same type) and pass `Opacity(typeof(T), n)` if the opacity
arrays are used. The coefficient and level data stay `Float32`/`Float64` and are widened as
needed. `test/forwarddiff_tests.jl` compares the derivatives with finite differences for every
ported type (except 4, which has no `rate` method). Notes:

- `to_single` and `to_double` (`constants.jl`) round the plain floating-point numbers to the
  precision of `ucalc`'s arithmetic and leave dual numbers alone; type 99, which mimics single
  precision, therefore uses `Float64` dual arithmetic when differentiated.
- Interpolation tables have kinks at their nodes (the derivative is that of one side), and the
  integrals of types 70 and 99 choose their step adaptively, so derivatives there are only as
  smooth as the rate itself: the finite differences of type 70 are noisy at the 1% level.

## Open items

- Type 85 (`PhotoionizeFeKedge`) sums sharp resonances, so it is sensitive to the
  last digit of the energy conversion: it uses `13.605692` as the single-precision
  literal that `ucalc` has (`Ry_eV_single`); with the double value the rates differ
  by 2e-5. `ucalc` also takes the resonance charge as the ion index minus 114 (an
  earlier ion numbering), which gives meaningless charges for the current ion
  numbers; reproduced as is.

- `ParPhotoIonize3` supports the 6-real form of the record, which is the only one in
  `atdb.fits`; `ucalc` also reads a 9-real form with two more parameters.
- `PhotoionizeDelta` stored its energies and line strengths one entry off before the real
  layout `[E_inf, E(1..m), f(1..m)]` was checked; fixed with the port of type 74.
- The collision rates (types 51, 56, 63, 98) return `fenergy` and `ienergy`, the forward and inverse
  rates times ΔE (erg cm⁻³ s⁻¹, `ucalc`'s `ans6` and `ans5`). Unlike the photoionization energies
  they are positive.

- For photoionization leaving an excited parent level, `ucalc` reads the energy of the
  final level (`rlev(1, idest2)` beyond the ion's levels) from stale memory; `fenergy2`
  and `ienergy2` use the continuum energy plus the parent level's energy instead and are
  only compared with `ucalc` for the ground parent level.
- 18 `ParPhotoIonize2` records (tables starting hundreds of eV above the shifted
  threshold) give negative photoionization rates; XSTAR gives the same values.

- `CollisionIonize` stores `ρ` from `rvec[2:end]`; the header says it starts at
  `r3`.
- The `Atom` records for iron and zinc look wrong (abundance 3.48 and 9.32, mass
  3.70 and 1.0); check the source table.
- `AtomicLine` (ucalc 4) has no `rate` method yet (an unfinished draft was removed); `atdb.fits` has no
  type 4 records to check a port against.
- The integer meanings of `ChargeExHp` are inferred.
- `DielecRecombH` and `TotDielecRecomb` store no data yet.
- `AtomicLine2`: the decay rate is `A` times the escape probability `pesc` and the
  photoexcitation rate (`irate`) comes from the `radiation` at the line energy, as in
  `ucalc`; the line opacity XSTAR adds to its continuum arrays (`linopac`) is not
  included.
- The second real of `ElectronImpact2` (`gf`) is not used by `ucalc`; its meaning is a guess.
- For type 98 with kind 1 or 4 and `C < 1`, `ucalc` can read its spline abscissa array
  out of bounds at low temperatures; no record in `atdb.fits` has that combination and
  Radix extrapolates the first segment.
- XSTAR's matrix assembly for bound-bound collision rates (lrtyp 3) decides which
  level is lower with `e1/e2 − 1 < 0.01` on the absolute level energies. For
  `ElectronCollision` (ucalc 51), where `idest1` is the upper level, this swaps
  the roles of nearly degenerate levels (energies within 1%). Radix always
  returns `init` = lower, `final` = upper, so the quirk is not reproduced.
- `RadiativeProb` (ucalc 54) computes the hydrogenic decay rate itself (Gordon's formula, `anl1`);
  the tabulated `A` of the record is zero and unused. Its energy output multiplies the rate by the
  dimensionless ΔE/kT instead of ΔE (`ener=delt` in the source), which looks like a bug; reproduced
  in `ienergy`. 60 of the 1,266 records link levels of equal `n` and give nothing.
- `EffectiveCharge` (ucalc 57) does not use the stored `Zeff`: it derives an effective charge from the
  level's ionization potential and returns the hydrogenic collisional ionization (`irc`, `szirc`) and
  three-body recombination of the level. 1,690 of the 5,608 records have five integers (no atomic
  number: `n, L, 2J, level, ion`); the constructor mis-assigned them before and now leaves `Z = 0`.
  `irc` has a separate formula for an effective charge of exactly 1, which no record reaches, so it
  is not checked against `ucalc`.
- `IronKAuger` (ucalc 86) returns the first Auger width (`A_widths[1]`, the second real of the record)
  as `frate` whatever the conditions; the radiative width and the width to the parent level are not
  used by `ucalc`.
- `RadiativeSuper` (ucalc 71) reads its table with the temperature running fastest (`A` has size
  `(nT, nₑ)`, stored as log₁₀ A; a one-point table holds A or log₁₀ A) and interpolates linearly in
  log₁₀ T (up to one dex outside the table) and log₁₀ nₑ. `ucalc` takes the statistical weight of the
  line from the first level of the record (`transition.lower`, which the header calls the lower
  level) as the upper one, and returns the decay as `irate`; reproduced as is. The decay of ions 96
  and 97 is capped at 10¹⁰ s⁻¹ and its energy uses a single-precision erg per eV. The line opacity
  `linopac` is not part of this type.
- `CollisionProb` (ucalc 63), transitions that change `n`: XSTAR's `ans1`/`ans2`
  put the large de-excitation-sized rate on the *upward* transition for any
  ordering of the two levels (the source comments "check if ans1 and ans2 are
  correct or inverted" and swaps them when `nf > ni` or `lf > li`). 3,597 of the
  6,015 records are affected. Radix reproduces ucalc literally (returns the
  levels in stored order, `frate` = `ans1`); worth reporting to the XSTAR
  maintainers.
- `ElectronCollision` clamps Υ at 0, which ucalc does not do: the spline gives
  small negative rates (about −1e-7) for 4 records.

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
| `ElectronImpact2` (98) | `transition`, `kind`, `ion`, `ΔE`, `gf`, `C`, `x_grid`, `Υ` |

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

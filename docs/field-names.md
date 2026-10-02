# Field names for Radix rate types

Applied in `src/`. The "Verified against the Fortran" section at the end records
what was checked and where the original proposal was corrected. The current names come from
the XSTAR data-type tables (`i1…i8`, `r1…`, `ionN`, `kN−1`), which are
positional labels and are inconsistent across types (`N` means ion number in
one type and electron number in another; the same level is `i`, `k`, `j`, `n`,
`iN`; `Nm`, `kNm`, `ionNm`, `NN`, `N_` all mean "the N−1 ion/level").

## Rules

1. **Name the meaning, not the Fortran slot.** Words for roles, single letters
   only for standard physics symbols (`Z`, `n`, `A`, `λ`, `σ`, `E`, `T`).
2. **One name per concept** across all types (table below).
3. **Rename the type code field.** `rate::Int8` shadowed the function `rate`.
   It holds XSTAR's *rate type* (`lrtyp`, e.g. 5, 8, 13) and is not redundant
   (`ChargeExH0` branches on it), so it is now `rtype`.
4. **Concrete, typed containers.** `Tuple`/`c::Tuple` are abstract and slow;
   use `Vector{R}` for tabulated data and fit coefficients (or `NTuple{N,R}`
   where N is fixed). `Tuple{R}` in `DielecRecomb2` is a one-element tuple and a
   bug. Integers are `I<:Integer`, floats `R<:AbstractFloat`.
5. **Units in a comment/docstring, not the name.** Energies eV unless a name
   ends in `_Ry`. Temperatures K. Cross sections Mb. Rates s⁻¹ or cm³ s⁻¹.
6. **Greek only where it is the standard symbol and easy to type**
   (`λ σ α γ ρ`). Avoid look-alikes: `Υ` → `upsilon`, `E∞` → `E_inf`.
7. **No trailing-underscore or mixed-case-suffix names** (`N_`, `Nm`, `S2p`).

## Shared vocabulary

| Concept | Name | Replaces |
|---|---|---|
| atomic number | `Z` | `Z` |
| ion index (XSTAR `ionN`) | `ion` | `N`, `ionN`, `N_`, `NN` |
| ion index of the N−1 (parent/target) ion | `parent_ion` | `Nm`, `ionNm`, `Nm` |
| single level index | `level` | `i`, `n` (AtomicLevel), `iN` |
| transition lower / upper level | `lower` / `upper` | `i`/`k`, `i`/`j` |
| level of the N−1 ion left behind | `parent_level` | `kNm`, `km`, `k` |
| principal quantum number | `n` | `n` |
| orbital angular momentum (term) | `L` | `L` |
| orbital quantum number of a subshell | `l` | `L` (ParPhotoIonize3) |
| spin multiplicity 2S+1 | `spin_mult` | `S2`, `S2p` |
| 2J | `twoJ` | `J2`, `J2p` |
| statistical weight 2J+1 | `g` | `J2p` (when it is 2J+1) |
| effective quantum number | `n_eff` | `ν` (comment wrongly says frequency) |
| transition type flag | `kind` | `type`, `t`, `it` |
| energy | `E` | `E`, `Ei`, `e` |
| threshold / ionization / limit energy | `E_th`, `E_ion`, `E_inf` | `Eth`, `Ei`, `Einf`, `E∞` |
| transition energy | `ΔE` | `ΔE` |
| wavelength | `λ` | `λ`, `e` |
| Einstein A (generic) | `A` | `A`, `a` |
| autoionization / radiative A | `A_auto`, `A_rad` | `Aa`, `Ar`, `cai` |
| weighted oscillator strength / oscillator strength | `gf` / `f` | `gf`, `f` |
| temperature | `T` | `T`, `t` |
| electron density | `ne` | `ne` |
| tabulated grids | `T_grid`, `ne_grid`, `E_grid` | `T`, `Te`, `ne`, `E` (tuples) |
| cross section | `σ` | `σ` |
| fit coefficients (unnamed) | `coeffs` | `c::Tuple`, `ρ`, `cfe` |
| named fit parameters | `a b c d e` (+ `T1 T2` break temps) | same |

## Per-type proposals

`label::String` is first in every type and omitted below.

**Levels and lines**

| Type (code) | Proposed fields |
|---|---|
| `AtomicLevel` (6) | `level, Z, ion, n, L, spin_mult, E, g, n_eff, E_inf` (`i`→`level`; drop duplicated "initial level" meaning) |
| `AtomicLevelFe` (83) | `level, ion, E, g` (`J2` was 2J+1 in the Fortran text; verify) |
| `AtomicLine` (4) | `lower, upper, λ, gf, A` (drop `e, f, a`; the in-progress `rate()` also needs level `E`/`g`) |
| `AtomicLine2` (50) | `lower, upper, Z, ion, λ, gf, A` |
| `RadiativeAPED` (91) | `lower, upper, Z, ion, λ, A` |
| `RadiativeProb` (54) | `lower, upper, Z, ion, A` |
| `RadiativeFeDecay` (82) | `lower, upper, ion, λ, E, gf, A_rad, A_auto` (`g` was really `gf`) |
| `TwoPhotonDecay` (76) | `lower, upper, ion, A` (`l` is the constant 1: drop) |
| `TwoPhotonRad` (5) | `lower, upper, a, b` (rename once the fit meaning is known) |
| `EffectiveCharge` (57) | `level, Z, ion, n, L, twoJ, Zeff` |

**Autoionization**

| Type (code) | Proposed fields |
|---|---|
| `Autoionize` (3) | `A, E` (`cai`→`A` prefactor, `eai`→`E`; units differ from a rate, say so) |
| `AutoionizeSat` (72) | `level, parent_level, Z, ion, spin_mult, L, g, A_auto, E` |
| `AutoionizeFe25Sat` (75) | `level, parent_level, ion, parent_ion, A_auto, E` (drop duplicated `N`/`N_`) |
| `IronKAuger` (86) | `parent_level, level, Z, parent_ion, ion, E, A_auto_total, A_auto_parent, A_rad` (currently `A::Tuple`; give the 3 values names) |

**Recombination and charge exchange**

| Type (code) | Proposed fields |
|---|---|
| `RadRecomb` (1) | `ion, A, η` |
| `TotRadRecomb` (38) | `Z, parent_ion, M, W, ion, A, B, T0, T1, C, T2` (`Nm`→`parent_ion`; keep Badnell's `M`, `W`) |
| `TotRadRecombH` (30) | `Z, ion` (needs `nmax`, currently undefined) |
| `DielecRecomb1` (7) | `ion, A, B, T0, T1` |
| `DielecRecomb2` (8) | `ion, C::Vector, E::Vector` |
| `DielecRecombH` (22) | `ion, a, b, c, d, e, f` (currently stores nothing) |
| `TotDielecRecomb` (39) | `Z, parent_ion, M, W, ion, C::Vector, T::Vector` (currently stores nothing) |
| `PhotoRecombX` (99) | see superlevel table |
| `ChargeExH0` (2) | `ion, a, b, c, d, T1, T2, ΔE` |
| `ChargeExHe` (9) | `ion, level, a, b, c, d, T1, T2, ΔE` (`i, k`→ recombining ion index and target level; confirm) |
| `ChargeExHp` (10) | `ion, level, a, b, c, d, e, T1, T2, ΔE` |

**Collisions**

| Type (code) | Proposed fields |
|---|---|
| `CollisionHlike1/2` (60/62) | `lower, upper, ion, coeffs` |
| `CollisionHelike` (68) | `lower, upper, Z, ion, coeffs` |
| `CollisionHeFine` (66) | `lower, upper, Z, ion, coeffs` |
| `CollisionHelikeSat` (73) | `lower, upper, Z, ion, coeffs` (constructor currently builds `CollisionHelike`; bug) |
| `CollisionLS` (69) | `lower, upper, Z, ion, coeffs` |
| `CollisionFe19` (81) | `lower, upper, Z, ion, upsilon` |
| `CollisionProb` (63) | `lower, upper, Z, ion` |
| `CollisionAPED` (92) | `lower, upper, Z, ion, T_grid, upsilon` |
| `CollisionIonize` (95) | `level, ion, E_th, T0, ρ::Vector` |
| `ElectronCollision` (51) | `kind, lower, upper, Z, ion, ΔE_Ry, C, upsilon_red::Vector` |
| `ElectronImpact1` (56) | `lower, upper, Z, ion, T_grid, upsilon` |
| `ElectronImpact2` (98) | `kind, lower, upper, ion, ΔE_Ry, C, upsilon_red::Vector` |

**Photoionization**

| Type (code) | Proposed fields |
|---|---|
| `ParPhotoIonize1` (49) | `n, L, twoJ, Z, parent_level, parent_ion, level, ion, E_grid, σ` |
| `ParPhotoIonize2` (53) | same as above |
| `ParPhotoIonize3` (59) | `ion_shell?` (`NN`, verify), `n, l, parent_level, parent_ion, level, ion, E_th, E0, σ0, ya, P, yw` |
| `PhotoionizeFeKedge` (85) | `n, L, twoJ, Z, parent_level, parent_ion, level, ion, Zeff, E_th_Ry, f, γ, scale` |
| `PhotoionizeDamp` (88) | `n, L, twoJ, Z, parent_level, level, ion, E_grid, σ` (the `k` and `i` fields in the current struct need checking) |
| `PhotoionizeDelta` (74) | `n, L, spin_mult, Z, parent_level, parent_ion, level, ion, E_inf, E_grid, f` |

**Superlevels (shared shape: a table over `ne_grid × T_grid`)**

| Type (code) | Proposed fields |
|---|---|
| `RadiativeSuper` (71) | `lower, upper, Z, ion, λ, ne_grid, T_grid, A::Matrix` |
| `CollisionSuper` (77) | `lower, upper, Z, ion, λ, ne_grid, T_grid, C::Matrix` (currently `A`) |
| `PhotoionizeSuper` (70) | `n, L, spin_mult, Z, parent_level, parent_ion, level, ion, ne_grid, T_grid, logα::Matrix, E_grid, σ` |
| `PhotoRecombX` (99) | same as `PhotoionizeSuper` with `α::Matrix` |

Drop the stored sizes `nd, nt, nx` (they equal the lengths) and store the 2-D
tables as matrices, not tuples of tuples.

**Atom / Ion**

| Type | Proposed fields |
|---|---|
| `Atom` (13) | `Z, n_electrons, mass, E_ion` (`N` here is electron count, not an ion index) |
| `Ion` (14) | `Z, ion, n_electrons, E_ion` (`I::I` field is the same name as the type parameter; rename) |

## Things to verify against the Fortran before renaming

- What `ionN` means numerically (ion index vs electron count) and whether
  `parent_ion` is always `ion − 1`; if so, drop that field.
- `AtomicLevelFe`/`AtomicLevel`: whether `J2p` is 2J or 2J+1.
- `ParPhotoIonize3` `NN`, `ChargeExHe`/`ChargeExHp` `i, k`, `PhotoionizeDamp` `k`.
- Whether `Float32` should be widened to `Float64` at load time (all physics is
  Float64 in XSTAR's double precision paths; the FITS data is single).

## Verified against the Fortran and the data

Checked against `ucalc()` in `xstarsub.f` and by loading `atdb.fits`:

- **`ion` is the global ion index** in the ion table (h_i=1, he_i=2, he_ii=3,
  li_i=4 … fe_xxvi=351), not an electron count. `Ion.N` was this index and
  `Ion.I` is the stage (1 = neutral), so they are now `ion` and `stage`.
- **`parent_ion` is `ion + 1`** for 287,173 of 287,180 `ParPhotoIonize1`
  records (the remainder are special cases), so it is kept as a field, not
  dropped.
- **`TotRadRecomb.Nm` is an electron count (N−1), not an ion index** →
  `parent_electrons`. `ParPhotoIonize3.NN` is the electron count of the ion →
  `n_electrons`.
- **`J2p` in `AtomicLevel`/`AtomicLevelFe` is 2J+1** (2P1/2 → 2, 2P3/2 → 4), i.e.
  the statistical weight → `g`. Other types store 2J as an integer → `twoJ`.
- **`rate::Int8` is not the data-type code** but the record's rate type
  (`lrtyp`), so it was renamed `rtype` and kept (correcting the first
  proposal, which said to drop it).
- **Charge exchange:** type 9 (`ChargeExHe`) takes `idat(1)` = own level and
  `idat(2)` → `nlevp + idat(2) − 1`, a level of the recombined ion; with a
  single integer both default to 1. Type 10 (`ChargeExHp`) ignores its integers
  in `ucalc`; in the data they are (level, ion index of the neutral, e.g. 29 =
  O I). Both are inferred, not documented.
- **`PhotoionizeDamp.k`** is a superlevel index present only in the 8-integer
  variant and 0 otherwise → `superlevel`.
- **`Atom`:** the two integers are both `Z` in every record (renamed `n_ions`,
  `Z`), and the first float is the abundance (H 1.0, He 0.1, C 3.7e-4), the
  second the atomic mass, not the other way round. The records for Fe and Zn
  (abundance 3.48, 9.32; mass 3.70, 1.0) look wrong; check the source table.
- **`TotRadRecombH`:** the Fortran uses the first integer (`Z`) as the effective
  charge, so `rate` uses `coef.Z^2` (my earlier `coef.N` guess was wrong).

Bugs found and fixed while renaming: `ratemap[9]` pointed at `ChargeExH0`
instead of `ChargeExHe` (59 records); `CollisionHelikeSat` records were built as
`CollisionHelike` (427 records); `TwoPhotonRad` called an undefined
`TwoPhoton`.

Left alone, but suspicious: `CollisionIonize` stores `ρ` from `rvec[2:end]`
(header says it starts at `r3`), `ChargeExH0.rate` clamps with `max(0, …)` and
tests `T > 5`, neither of which is in `ucalc`.

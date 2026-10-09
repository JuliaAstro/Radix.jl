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
   `ntot` (`xpx`, the hydrogen density). The level data and element masses are stored in the
   coefficients (see Level data in the coefficients), and the radiation field is a separate
   object.
6. **Greek only for standard, easy-to-type symbols** (`λ σ α γ ρ η Υ`). Only
   `E∞` is spelled out (`E_inf`), because a subscript is not a valid identifier
   character.
7. **Group values that travel together in small structs** (below).
8. **No numeric literals in expressions.** Every empirical or physical number is a
   named constant. Physical constants come from PhysicalConstants.jl through
   `Radix.constants()` (see Physical constants); other numbers used by more than one file are
   in `src/constants.jl` (`cx_unit`, `tiny`, `Mb`, …);
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
  integer both default to 1. Type 10 (`ChargeExHp`) takes its first integer as the level
  (for the record with one integer that is the ion index, e.g. 29 = O I); the others are
  (level, ion index of the neutral). Both are inferred.
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
spectrum on it (`bremsa`). The rates take it as `radiation=` (`NO_RADIATION` by default: no photons, on
XSTAR's standard grid `xstar_energy_grid()`), the escape
probabilities as `ptmp=` and the populations as `abund=`. `Opacity(n)` holds the
continuum opacity and recombination emissivity arrays that the integrals add to when
passed as `opacity=`. The shared engine is `photoionize_level` (types 49 and 53) with
`photoionization_integrals` (XSTAR's `phint53`) and `phextrap`.

## Level data in the coefficients

A rate needs more than its own coefficients: the energies, weights and quantum numbers of the levels it connects,
the number of levels of the ion (the continuum is the last) and the atomic mass of the element. These come from a
`Levels`, which `levels(records)` builds from the `AtomicLevel` and `AtomicLevelFe` records (the Fe UTA
levels, which hold only an energy and a weight, are added as levels with zero quantum numbers) and the `Atom` and
`Ion` records (the mass of each ion's element). It behaves as the dictionary `(ion, level) => AtomicLevel` and also
answers `nlevels(table, ion)`, `atomic_mass(table, ion)` and `level_counts(table)`.

The rates that need it have a last field `levels` (a `Levels`), and the constructor of their database record takes the
table as an extra argument. `load` builds the records of the levels, elements and ions first, makes the table from them
and then constructs all the other records with it (`construct(T, rate, label, ivec, rvec, levels)` passes the table to the
types for which `needs_levels(T)` is true and ignores it for the rest), so `rate(coef, cell)` is all that is needed. A coefficient
made by hand takes the table as its last argument. The number of levels and the mass (`AtomicLine2`, `RadiativeFeDecay`,
`RadiativeSuper`; `mass=` still overrides it) are derived from the table, a table without the levels of the ion gives no
rate, and no keyword of a rate is left without a default. The photoionization rates default to `radiation=NO_RADIATION` (no
photons).
Exceptions: `ChargeExH0`, `ChargeExHe` and `ChargeExHp` keep their `nlev=0` keyword, because their records do not identify
the ion reliably (`ChargeExHe` has no `ion`, and `ChargeExHp` has one record with a single integer).

The `ucalc` tests build their table from the levels of the fixture and the atomic mass of the driver, and
`test/levels_tests.jl` checks the table, the construction and the defaults.

## Physical constants

The rates take their physical constants from `Radix.constants()`, a `Constants` object evaluated
from PhysicalConstants.jl with CODATA 2022 by default: k, h, c, the Rydberg, erg per eV, hc/k and
its temperature floor ΔE/50k, Ry/k, eV/k, the proton-electron mass ratio, the line cross section
π r_e c, the thermal speed of 1 amu at 10⁴ K, and the Maxwellian collision and Saha prefactors. Use
another CODATA set with `set_constants!(Constants(PhysicalConstants.CODATA2018))` (also CODATA 2014),
or temporarily with `with_constants(f, K)`; `Constants(K; field=value)` copies a set with changes.

XSTAR's `ucalc` rounds these differently, sometimes in different ways within one program
(`kT_eV` = 0.861707, 8.617e-5 and 1.38066e-16 erg/K for k, three values of hc/k, 1800 for the
proton mass). `ucalc_constants()` has its numbers, field by field; the tests that compare with
`ucalc` (`ucalc_tests.jl` and the formulas transcribed in `runtests.jl`) run with that set, and
`constants_tests.jl` checks that the CODATA values are consistent with each other and within the
rounding of `ucalc`'s (1e-3, 3% for the proton mass). Where `ucalc` is sensitive to the rounding (the
e^{-ΔE/kT} factors) the CODATA rates differ from it by more than the 1e-6 of the tests: the
Boltzmann factor amplifies the 3e-5 difference in k by ΔE/kT.

4π and 8π (the normalization of the spectrum in the photoionization integrals) and π (the resonance
profiles of type 85) are the exact values, and `ucalc` has 12.56, 25.3 and 3.14159. The photoionization
rate is independent of the choice, since the factor multiplies and divides the flux, but the
recombination and emissivity terms are not: `ucalc`'s 25.3 is 0.7% above 8π.

Coefficients that XSTAR writes as a number but that are combinations of constants are derived too, and
checked against its number (within 1.3e-3, except the Saha factor of type 57, 3.5e-3):
`bb_coeff` = 1/(h³c²) and `fo_saha` = Saha prefactor × 8π (the photon density and the Saha factor of
the recombination integrals), `milne_coeff` = 4π/((2π m_e)^{3/2} c²) × 1 Mb, `A_to_f` = 1/(8π² r_e c),
`gordon_A` = (2π/3) α³ c R∞ with the reduced mass of hydrogen (with it the 2p → 1s rate is 6.2649e8 s⁻¹,
as NIST has it; the infinite-mass value is 5e-4 higher), `ps_alfa` = m_e/2k and `ps_pd` the Debye length
coefficient of the l-changing collisions, `Ry_K_sz` = Ry/k and `sz_rate_coeff` of Simpson and Zhang's rate,
and the Saha factor of type 57. `constants_tests.jl` also checks the two routes to the recombination
coefficient, the Milne integral and the photoionization code's, against each other (1-2%, the accuracy of
the integrals).

Not derived: XSTAR's fit coefficients and numerical conventions (the fits of `impcfn`, the Bethe and
Simpson-Zhang tables, `ps_rho_coeff`, `ps_b_coeff`, `impactn_psi_coeff`, `impactn_cr_coeff`, `szirc_coeff`,
`erc_f_coeff`, the ionization fit constants of type 57, `pexs_a_coeff`, ...). `irc_coeff` = `erc_s_coeff` =
1.095e-10 is within 1.8e-3 of 2 (8k/πm)^{1/2} π a₀² per √K, which is probably its origin, but this is
unconfirmed.

## Heating, cooling and thermal equilibrium

`src/heating.jl` and `src/mixture.jl` (`element_heating`, `heating_cooling`, `thermal_equilibrium`) port XSTAR's `msolvelud`/
`msolvelucy` (the sums over the populations), `calc_hmc_all`, `comp2`/`cmpfnc`, `freef`, `bremem` and `heatf`, and the
iteration of `dsec`. What the code does and what was found:

- **Two channels.** ucalc returns two pairs of energies: `ans3`, `ans4` (the energy of the photons, "photon pov": the heating
  is what the radiation field deposits and the cooling what the gas radiates, so that its imbalance is what the temperature is
  iterated on) and `ans5`, `ans6` (the energy of the electrons). Each record puts `ans4` at its lower level and `-ans3` at its upper
  one (`ans6` and `-ans5` in the second channel); a positive entry times the population of its level is a cooling, a negative one
  a heating. `ucalc_energies` gives the four from the energies of the result of `rate`, by type, and the signs are those of ucalc:
  the lines and `PhotoionizeSuper`/photoionization negate them (`ans3 = -ienergy`, `ans4 = -fenergy`, `fenergy2`, `ienergy2` in the
  second channel), `RadiativeFeDecay` (type 82) does not (`ans3 = ienergy`, `ans4 = fenergy`), the radiative decays 54, 71 and 76
  have only `ans3`, and the collisions only `ans5` and `ans6`. Collisional excitation cools through the line that radiates it.
- **`hmctot = 2 (ht - cl)/(1e-37 + ht + cl)`** and the temperature is iterated until `|hmctot| < 1e-4` (steps of 1.2, twice that if
  `|hmctot| > 0.9`, then false position with the unmoved end halved). The imbalance changes by only 0.0035 per 10⁴ K near
  the equilibrium of the reference run, so a relative error of 0.1% of heating or cooling is a 2% error of the temperature.
- **`bremsmap` and the flat tail** (`map_spectrum`). XSTAR does not use the spectrum on its grid but `bremsam(m) = bremsa(nbinc(epim(m), epi))`,
  and `nbinc` never returns a bin above `n - max(2, n/50)`. With the same grid every bin above it takes the flux of that bin: a flat
  tail where a power law falls. The Compton heating, which the high energies dominate, is 50% larger than for the unmapped spectrum for an E⁻¹
  spectrum on the grid of 999 bins. Radix does not map the spectrum itself: pass `map_spectrum(radiation)` to reproduce XSTAR.
- **Compton.** `coheat.dat` (a table of 101 × 101 values, read by `load(Compton, path)`) is in the data directory of XSTAR, next to `atdb.fits`.
  `integral(Compton, radiation, T)` (the `comp2` integral; `σ(compton, Eph, Te)` is `cmpfnc`) is checked against the real F90 `comp2` (`test/reference/ucalc/drvcomp.f90`) to 10⁻⁶. The processes that
  are not tied to levels are types (`src/processes.jl`, subtypes of `AbstractContinuum`, whose coefficients are the functions `α` (absorption), `j` (emission) and `σ` (scattering)) with generic functions: `heating(process, radiation, T, nₑ)`,
  `cooling(...)`, `heating_cooling(...)` (both) and `α(process, E, T)`, `σ(process, E, T)`, the absorption or scattering coefficient, since each coefficient is a function of the process (the opacity is that times the density, `opacity(process, E, T, nₑ)`); `FreeFree()` (absorption: heating and opacity),
  `Bremsstrahlung()` (emission: cooling), `Compton(Te, Eph, σ)` (heating and cooling; the energy exchanged per electron is the coefficients `heating(compton, E, T)` and `cooling(compton, E, T)`, from which the heating and cooling derive) and `Thomson(cfrac)` (an opacity, which is also in the
  continuum opacity `opakcont`: `AbstractScattering`). A process without a heating, cooling or opacity has 0. `standard_processes(compton)` makes the four
  of XSTAR; `heating_cooling(mixture, ...)`, `thermal_equilibrium` and `march_zones` take such a collection and report the term of each
  (`hc.processes`). The free-free terms use XSTAR's fitted coefficients (Gaunt factor 1, ion density 1.4 nₑ).
- **Photoionization that leaves an excited level of the next ion**: ucalc reads the energy of the final level from stale memory, and
  the energy measured from the threshold is a difference of nearly equal terms (it reached 10⁵ times the cooling). These records are left
  out of the second channel (the first, from zero, is verified for every parent). XSTAR keeps finite values for them, so the second channel
  of the heavy elements (S, Ar, Fe, ...) is 10–90% low.
- **Which database.** The reference outputs were made with XSTAR 2.59j, whose database (`$HEADAS/refdata/atdb.fits`, 874.7 MB) is not the
  one in the source tree (`ftools/xstar/data/atdb.fits`, 871.1 MB). For hydrogen the older one has records of type 60 (the 2s-3s record has
  the levels 1 and 7, which makes 3s 7.5 times too populated and the cooling of hydrogen 6-10% too large) and the newer has none. With the
  database of the package every element agrees with XSTAR's table of the first zone to 0.0-0.2% in both channels (the same bias in heating and cooling,
  which comes from XSTAR's Lucy iteration, converged to 1% only) and the electron fraction and the temperature of the equilibrium
  to 10⁻⁵ and 2×10⁻⁴; with the database of the tree the temperature is 4-5% low. The Compton and bremsstrahlung terms agree to 5×10⁻⁴.
- **Records from the ground level to a level of the next ion.** The records whose rate type is ordered by the energies of the two levels (all but 7 and 41) are put lower-first by `calc_hmc_ion`
  with the energies of `leveltemp`, which for a level beyond the ion (`final > nlev`) are stale ones from the ion before. Radix left such records out (`l = u = 0`), but one from the ground level is never switched (its energy is 0), and the
  `ParPhotoIonize3` of rate type 1 of Fe VII and Fe VIII (to the levels 3 and 4 of the next ion; 5% of the photoionization of the ground level of Fe VIII at log ξ = 1) are in the matrix. Without them the iron
  of the first zone of the thick slab with the iteration of the temperature was 4-7% off in the ions (heating of iron -7%, cooling -2.4%) and T 2.9% low (8% with the database of the tree); with them the ions agree to 10⁻⁵, iron to 4×10⁻⁵
  and T of all the zones to 0.05% (`test/reference/xstar_fe_balance`). A record from an excited level to a level of the next ion is still left out.
- **Not included**: the escape probabilities from the optical depths (the second zone of the reference run has line trapping:
  its equilibrium temperature is 0.5% low), and the transfer.

## Transfer: escape probabilities and optical depths

`src/transfer.jl` ports `pescl`/`pescv`, the opacities of `calc_emisab_ion`, `stpcut` (the accumulation of the depths) and the use of
the depths in `calc_hmc_ion`:

- **What is accumulated.** Each line (rate types 4 and 9) has a centre optical depth `tau0` and each recombination edge (rate type 7) a
  `tauc`, in two directions; `stpcut` adds `opacity × Δr` after every zone (`add_zone!`). The opacity of a line is the cross section at
  its centre (the `opacity` of `rate`, ucalc's `opakab` without the abundance) times the density `ntot × abundance × population` of its
  lower level; that of an edge is the `opacity` of the photoionization record called with `abund = (population × abundance of its
  initial level, of its final level)`. The next zone's escape probabilities follow: `(pescl(τ₁)(1-c), pescl(τ₂)(1-c) + 2 pescl(τ₁+τ₂) c)`
  for a line (`c` the covering fraction), `pescv` for an edge, 1/2 each for the others (`escape_probabilities`). The continuum
  depth `dpthc` (the attenuation `exp(-dpthc)` of the spectrum in `trnfrc`) is in the next section.
- **Single-precision literals in `pescl`/`pescv`**: 1.2, 1.e-5 and 1.e-12 are single precision, which makes the double-precision
  values differ by 4×10⁻⁸; checked against the real functions (`test/reference/ucalc/drvpesc.f90`) to 10⁻¹³.
- **`trnfrc` writes the flux as `L/(12.56 r²)`, not `L/(4π r²)`**: 0.05% smaller. `point_source` uses the `fourpi` of the constants
  in use (4π by default, 12.56 with `ucalc_constants()`), which is what the 0.1% bias between Radix and XSTAR's heating and
  cooling was.
- **The reference run has two zones.** `nsteps=2` is a slab of two zones of 5×10¹² cm: the first computed at 10¹³ cm (`frac_heat_error`
  and the first row of the table of the second), the second at 1.5×10¹³ cm (the third row), and the optical depths in `xout_lines1.fits`
  (`depth_inward`) and `xout_rrc1.fits` (its `depth_outward`, which is the inward one) are their sum. This reproduces the line depths of
  all 451 lines above 10⁻⁹ to 1% (median 0.9999984, with XSTAR's constants) and the edge depths of the H-like and He-like ions
  (the others have several records with the same level and threshold and no parent level in the table to tell them apart) to 0.1-0.3%.
- **Continuum opacity and attenuation** (`src/continuum.jl`, `Continuum`, used by `march_zones`). `calc_emis_all` starts
  `opakc` at Thomson scattering `nₑ σ_T (1 - cfrac)` (`Thomson(cfrac)` (the coefficient `σ_T (1 - cfrac)` times `nₑ`)), adds the photoionization opacity of the records (called with the populations as
  `abund`) and then free-free (`FreeFree()` (its coefficient times `nₑ nᵢ`), only in the total, not in `opakcont`); `stpcut` adds `opakc Δr` to `dpthc`. Only the
  `nrank = 10` strongest edges of each energy bin enter (`rlbin`, ranked by the edge opacity of the zone; the bin of an edge is that
  of its edge energy, `E_grid[1] × 13.598` for type 49 and `max(0.1, E∞ - E)` of the level for the others, set once by `xstarsetup`)
  and every record of rate type 42. The transmitted spectrum of the output (`writespectra3`) is `incident × exp(-dpthcont)`, with `dpthcont` the depth of `opakcont` (Thomson and photoionization: no lines or free-free), while
  the spectrum that passes the zones loses the `opakc` (lines and free-free included) of each (`transmit`; `march_zones` returns `dpthc` and `dpthcont` as well). Against XSTAR's slab of 10²¹ cm⁻²
  (`test/reference/xstar_thick_slab`) the depth `dpthcont` of all 999 bins agrees to a median of 1.0004 and within 5%
  (Thomson scattering 0.1%, the He II edge 2.5%).
- **Lines in the bins** (`src/lines.jl`, `add_line!`, `voigt`). `linopac` puts the opacity of a line (its centre opacity `oplin` times the Voigt profile of the thermal and turbulent widths
  and the damping parameter `delea/(width 4π)`, averaged in the bins of the grid) into `opakc`. `calc_emis_ion` calls it for the `nrank` lines of data type 4 with the largest
  emissivity `rcem` in each bin (`rlbin` with `lopak = 0`; `line_emissivities`, the energy of the decay split by the escape probabilities as `calc_emisab_ion` makes it); only the type 50 lines
  (`AtomicLine2`) call it. `rlbin` never stores an entry whose place is the last of the list (`rank!`, which the ranking of the edges uses too). Checked against the real `linopac` and
  `voigte` (`test/reference/ucalc/drvlinopac.f90`): `voigt` is identical and the bins of four lines agree to 10⁻⁷ (the single-precision `12.9` of the thermal speed), including the quirk that the profile
  loop never ends early (`ml1min` stays above the bin of the line, so the wings fill the grid out to 10⁴ steps). `delea` is `A × 4.136×10⁻¹⁵` eV; the records of type 41 of the database, which
  `deleafnd` looks for first, are not read. Not ported: the Fe UTA lines (type 82), whose `delea` is its `A_auto`.
- **Two-photon continua** (`add_two_photon!`). The decays of type 9 (`TwoPhotonDecay`, 85 records, and the four `AtomicLine2` of rate type 9) are not lines: `ucalc` spreads the energy `A hν`
  of the decay over the bins below `hν` with the shape `E² (hν - E)` (normalized by its integral from 0, which starts at the second bin) in `rccemis`, with the population of the upper level and `ptmp` of the optically thin line (`Continuum` adds them).
  Without them the diffuse emission of the reference slab misses the 2γ continua of the H-like and He-like ions that fill the bins between the edges (a factor 6 at 588 eV).
- **The spectrum that passes through the zones** (`Transmitted`, `transmit`: XSTAR's `heatt` and `trnfrn`). In the first pass (`ldir = -1`) `trnfrc` makes the radiation of each zone from `zrems(1)/(4π r²)`, and `heatt` updates it:
  `zrems(1) = max(0, zrems(1) - (F κ - 4π(ε₁+ε₂)) fac Δr 4π r²)` with `κ = opakc`, `fac = (1 - e^{-κΔr})/(κΔr)` (1 below 0.01), `ε₁ = rccemis₁ + brcems (1-c)/2`, `ε₂ = rccemis₂ + brcems (1+c)/2`;
  `zrems(2)` to `zrems(5)` accumulate `4π ε fac Δr 4π r²` (with `opakcont` for the last two), and the last two are the columns `emit_inward` and `emit_outward` of `xout_cont1.fits`. It is
  `exp(-κΔr)` for the attenuation (to first order in a thin zone), plus the emission, so that the field of the next zone is not `zremsz exp(-dpthc)`. **The bremsstrahlung of the source code is 4π too large** (it adds `brcems`, in all directions, to the emission per steradian
  of the recombination continua, and multiplies the sum by 12.56): the XSTAR 2.59j that made the references does not do it (bin 1 of `xstar_pow_xi2` agrees to 10⁻⁶ and the bremsstrahlung bins of the thick slab to 4×10⁻⁶ without the factor,
  and are 12.56 times too high with it), so `transmit` divides it by 4π. With the XSTAR zones of the thick slab (`delta_r` of its table: five zones of 1.8×10¹⁶ to 2.8×10¹⁶ cm)
  the emission columns of 571 bins agree to a median of 1.0002 and 1.0015 (inward and outward) and 5% at worst, and the He II and O VIII fractions at the inner edges of the zones to 2×10⁻⁴-9×10⁻⁴ and 3×10⁻⁵-1×10⁻⁴ (0.7-2.7% and 0.07-0.3% without the emission). The last row of the table (the end of the slab) is not a zone:
  XSTAR's final calculation uses the field of the last zone before `heatt`, which Radix does not do. Not ported: the luminosities of the lines and edges (`elum`, `elumab`), which `heatt` also updates, and the further passes (`ldir > 0`).
- **The thickness of the zones** (`march_slab`, `step_thickness`: XSTAR's `step`). `xstar` solves zone after zone from the distance `r` of the source until the column `xpxcol` is reached; the first zone has no
  thickness (it only gives `opakc`), and the thickness of the others is `min(xpxcol/xpx, r/numrec, (xpxcol - xcol)/xpx, emult/opakc)`, the last over the bins above `ectt = 1` eV whose depth `dpthc` before the zone does not exceed `taumax`
  and whose spectrum `zrems(1)` is above 10⁻¹² (10³⁸ erg s⁻¹ erg⁻¹), with the `opakc` of the zone before (lines included: the bin of the He II line at 41 eV limits the zones of the reference slab). `numrec` is `nsteps` (not a number of zones) and `emult` `emult`
  (the defaults are 0.75, `taumax = 5`). For the thick slab of 10²¹ cm⁻² (`emult = 1`, `nsteps = 3`) it reproduces the seven zones of XSTAR's log (log N = 20.26, 20.56, 20.74, 20.86, 20.95, 21.00 at their ends): their depths agree to
  10⁻⁴ (1.8321, 3.6643, 5.4553 and 7.2120 ×10¹⁶ cm against the 1.8322, 3.6644, 5.4550 and 7.2111 of the table of the run; the table has no row for the sixth zone). The final calculation at the end of the slab and the conditions on the temperature and the electron fraction of the loop of XSTAR are not done.
- **The luminosities of the lines and the edges** (`Luminosities`, `add_luminosities!`, `line_emissivities`, `edge_emissivities`: `elum` and `elumab` of `heatt`). Each zone adds `rcem Δr 4π r²` to the inward and outward luminosity of a line (the two escape
  shares `ptmp₁`, `ptmp₂` of its decay energy) and half of `(cemab₁ + cemab₂) Δr 4π r²` to each of an edge, in 10³⁸ erg/s; they are the `emit_inward` and `emit_outward` of `xout_lines1.fits` and `xout_rrc1.fits`. `cemab` is made for all the edges (the records of rate
  type 7 with a level of the ion above 10⁻³⁴), not only for the ones that `rlbin` ranks. For the thick slab with the zones of XSTAR the total of the edges agrees to 3×10⁻⁵ (807126 against 807151) and the 599 lines above 100 in its table to a median of 1.0001 and 5% at worst, but for
  hydrogen: its lines are 2× (Lyα) and 20× (the 2p-3s of Hα) off because the database of the tree has other records for it than the package of XSTAR that made the table (`xout_lines1.fits` of the reference is from the package).
- **The source and the model** (`src/spectrum.jl`, `power_law`, `source_distance`; `slab_model`). `ispec4` makes the spectrum `E^α` (`trad = α`; 10⁻²⁴ below 0.01 eV) and normalizes it to the luminosity of 1 to 1000 Ry in the bins
  `nbinc(13.6)` to `nbinc(1.36e4)`; `ispecgg` normalizes it again in the bins whose energies are in the range (with 13.6 as a single-precision literal), which is not the same bins. The result is the incident
  spectrum of the reference runs to 5×10⁻⁶ (the precision of the files). The distance of the source follows from the ionization parameter, `r² = L/(ξ n)`. `slab_model(mixture, processes; density, column, logξ, luminosity, α, ...)` puts it all together:
  `power_law`, `source_distance` and `march_slab`, returning the zones, the spectra and the `transmitted` spectrum `incident × exp(-dpthcont)` of the output. The turbulent speed `vturbi` only gives the widths of the lines
  (XSTAR's `gsmooth` of the continuum for `vturbi > 0` is not done: see the next item). The other spectra are `blackbody` (`starf`: `trad` in 10⁷ K, `1/kT = 1.16×10⁻³/trad` eV⁻¹, called by `xstar` with the luminosity 1), `thermal_bremsstrahlung` (`ispec`, `kT = 861.707 trad` eV) and `tabulated` (`ispecg`: a table interpolated in the logarithms, ascending in the energy, with the units of `spectun`); the first two agree with XSTAR's incident spectra to 3×10⁻⁵ (`test/reference/xstar_sources`).
- **Constant pressure** (`ConstantPressure`, `gas_density`, `source_distance_pressure`; `slab_model(...; pressure)`). With `lcpres=1` the hydrogen density is `xpx = p/1.38e-12/t` (T in 10⁴ K, the pressure being that of the hydrogen) and `calc_hmc_all`,
  `calc_emis_all` and `calc_emisab_all` recompute it for the temperature they are given, so the density follows the temperature of the iteration of each zone; `rlogxi` is then the ionization parameter of the pressure, `r² = L/(4π c P Ξ)`
  (with `c = 2.99792458e10` as a single-precision literal and 12.56 for 4π). Radix takes `ntot` as a number or a function of T (`thermal_equilibrium` evaluates it at each temperature), and the column of the zones adds up with the density that each has at its
  final temperature. The reference run (`test/reference/xstar_const_pressure`, `P = 10⁻⁷`) gives r to 10⁻⁵ and the first zone (T, n, x_e) to 4×10⁻⁴, 3×10⁻⁴ and 10⁻⁵ with the database of the package. **XSTAR does not always converge the temperature**: in its
  second zone it stays at the T of the first with a heating imbalance of 0.15%, which Radix reproduces (0.154% in the same state: the physics agrees), where the iteration of `dsec` stops when the false position repeats a temperature (`testt < 2×10⁻⁹`); Radix finds the
  root there (10.71 against 10.46, 2.4%), as it does for the zone 3 (+0.4%) of the thick slab with constant density. The other `lcpres` values are not ported.
- **A density that varies as a power of the distance** (`PowerLawDensity`, `slab_model(...; radexp)`; `lcpres=0`). The density of a zone is `xpx = xpx0 (r/r0)^radexp` at the radius of its inner edge, `r0` the radius of the first zone, but `xstar` updates it for the next zone
  before it adds the column of the zone, so **the column of a zone is `Δr` times the density at its outer edge**. Radix does this (the thickness of the next zone, `step`, uses the density at its position). For `radexp = -3`, `rlrad38 = 1e-2`
  the rows of the table of XSTAR (log r = 15.50, 15.50, 15.62; log N = 18.65, 18.84; log n = 4.00, 3.63) are reproduced (`test/transfer_tests.jl`). Such a slab has a finite column (n r⁻³ from 3×10¹⁵ cm adds to 1.6×10¹⁹ cm⁻²), and when it is asked for more, XSTAR goes on
  to 10⁵² cm until its buffer of 3999 steps is full; Radix, with the density going to zero (zones of a few cm in the end), has stopped with a singular matrix in a zone near the saturation (column 1.1×10¹⁹ cm⁻²), not looked into.
- **`gsmooth2` zeroes the continuum when `vturbi > 0`.** XSTAR 2.59j smooths `opakc`, `rccemis` and `brcems` with a Gaussian of width
  `E vtherm/c` whenever `vturbi > 1e-34` (the default is 1 km/s). That is about 10⁻⁴ E, and the bins of a grid of 999 points are 1.6% wide
  (of 9999 points 0.16%: the same for T below 4×10⁵ K): the loop ends at the first neighbour with `exp(-earg) = 0`, and never adds the bin
  itself, so the smoothed arrays are 0 below 20 keV. The continuum is then not absorbed in the slab and has no diffuse emission (the
  `transmitted` column of the output is the incident spectrum and `emit_inward` is 10⁻¹⁷ instead of 10¹⁵). Radix does not do it; use
  `vturbi=0` for reference runs of the continuum (the lines are unaffected).

- **Line trapping is negligible in that slab**: the largest depth is 4×10⁻³ (O VIII Lyα), which changes the heating and cooling of the second
  zone by 10⁻⁷. The 0.5% difference of the equilibrium temperature of that zone is not trapping.

## Reference tests

`test/ucalc_tests.jl` compares each ported rate with XSTAR's real `ucalc` on records
sampled from `atdb.fits` (`test/reference/ucalc/`, with the driver, the build script
and the generator). The driver is built like HEASoft's `xstar`, so XSTAR's
single-precision literals limit the agreement to about 1e-6. Types covered: 1, 2, 7, 8, 9, 10, 22, 30, 38, 39, 49, 50, 51 (5- and 9-point fits), 53, 54, 56, 57, 59, 60, 62, 63, 66, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 81, 82, 85, 86, 88, 91, 92, 95, 98 and 99.

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
  last digit of the energy conversion: with `ucalc_constants()` it uses `13.605692` as the
  single-precision literal that `ucalc` has (`Ry_eV_single`); with the double value the rates
  differ by 2e-5. `ucalc` also takes the resonance charge as the ion index minus 114 (an
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
- `DielecRecombH` and `TotDielecRecomb` store no data yet.
- `AtomicLine2`: the decay rate is `A` times the escape probability `pesc` and the
  photoexcitation rate (`irate`) comes from the `radiation` at the line energy, as in
  `ucalc`; the line opacity XSTAR adds to its continuum arrays (`linopac`) is `add_line!`,
  which `Continuum` calls for the lines that `rlbin` selects.
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
  `linopac` is not called for this type (XSTAR does not for it either).
- `CollisionSuper` (ucalc 77) reads its table like `RadiativeSuper` (temperature fastest, `C` has size
  `(nT, nₑ)` and holds log₁₀ of the de-excitation rate) through the shared `interpolate_super_table`. The
  excitation rate uses the wavelength of the record in the Boltzmann factor (`hc_over_k`; `ucalc` has
  1.43817e8) and a weight `2(2l + 1)` where `l` comes from numbering the first level of the record by
  shells (`level` 1 is 1s, 2 and 3 are 2s and 2p, ...), as `calt77` does with the level index. The
  temperature floor ΔE/50k, from the energies of the levels, is `T_floor_coeff` over the wavelength
  from the single-precision 12398.4016 of `ucalc` (`hc_eVÅ_single`): with `ucalc_constants()` the
  factor e^{-ΔE/kT} amplifies its 2e-8 rounding to 1e-6.
- `RadiativeAPED` (ucalc 91) is `ucalc`'s type 50 (it jumps to that code): `rate` builds the equivalent
  `AtomicLine2` and calls it, with the same keywords and results (the stored second real, 0, is not used).
- `TwoPhotonDecay` (ucalc 76) returns the stored `A` as `irate` (independent of the cell), from the higher to the lower
  level by energy, with `ienergy = A ΔE`. `ucalc` also adds the two-photon continuum (`E²(E_max − E)` normalized to `A`)
  to its emissivity array; that is not included, as for the line opacity of types 50 and 82.
- `ChargeExHp` (ucalc 10, four records) multiplies the H⁺ density, which `Cell` does not carry (it has the neutral
  hydrogen), so `rate` uses `ntot - nₕ` for it (`xh1`; the `ucalc` test builds its cell with `nₕ = xpx - xh1`). The reals are `a, b, c, d`, the range of the fit
  `T_min` and `T_max` (K, not applied) and `E_k` = ΔE/k in 10⁴ K (the old fields `e`, `T1`, `T2`, `ΔE` had these
  shifted by one); the eighth is not used. `ucalc` takes the level from the first integer of the record even when
  it is the only one (the ion, in one record) and does not look at the levels; `1 + c e^{dT}` is not limited at 0,
  unlike in type 2, and there is no temperature limit.
- `AutoionizeSat` (ucalc 72) and `AutoionizeFe25Sat` (75) use `calt72`, `3.3e-11 (13.6 eV/kT)^{3/2} e^{-E/kT} (A/10¹³) g`
  (`g` is 1 for type 75, whose records have two reals), a capture-like rate built from the stored autoionization rate.
  `ucalc` multiplies it by the electron density and calls the result the rate of the transition (`irate`); type 72 also
  returns `frate` as that times `nₑ`, the Saha factor `saha_ci g(1)/g(continuum)/T^{3/2}` and `e^{E/T}` with the energy in
  eV divided by the temperature in K (≈ 1: the units look wrong in the source), and type 75 computes the same product
  and discards it (`frate` = 0). These are as `ucalc` has them; the physics of the two rates is not checked. The 13.6 eV is
  `Ry_eV_coarse`, now also used for the thresholds of types 70 and 99, the ionization potential of type 57 and the
  extrapolation of cross sections (`ucalc` writes 13.6 in these places and 13.605692 in others).
- `CollisionLS` (ucalc 69), `CollisionHelike` (68) and `CollisionHeFine` (66) are three fits of the He-like effective
  collision strength with the same rate: type 69 (Kato and Nakazaki; 6 coefficients, or 9 with a non-resonant part,
  the first being an energy in eV that scales the temperature) is clamped at 0 and has no temperature floor; type 68 is
  a quadratic in log₁₀(T/Z³), clamped at 0, with the floor ΔE/50k; type 66 (fine structure; six coefficients in the
  database, `calt66` also reads two more terms) uses the energy of the record, not of the levels, for the Boltzmann
  factor, the floor and the energy of the rates, and is not clamped. `ucalc` reads the first coefficient of types 66
  and 69 in eV with `eboltz = 1.160443e4` K/eV (`eV_K_ls`). The scaled energy is limited to (0.05, 77) in type 69 and
  to 77 in type 66; at 77 the fits cancel two large terms (y²/2 times `d`), which amplifies any rounding of the
  exponential integral by ~3000: `ucalc`'s `expint` has single-precision coefficients (`single_expint` of
  `ucalc_constants()`, 2e-7 effect elsewhere), and the 66 and 69 fits differ from `ucalc` by 1e-4 without it. The
  default uses the full coefficients. `calt69` prints an error and leaves Υ unset for y < 1e-20, which the
  database does not reach.
- `CollisionHlike1` (ucalc 60) and `CollisionHlike2` (62, which `ucalc` runs through the same code) evaluate the
  fits of `calt6062` in the scaled temperature τ = kT/Ry (the constant `Ry_per_K`; `ucalc` has 6.33652e-6, which is
  Ry/k for the Simpson-Zhang value 1.578203e5 K), limited to 1 (above 10⁹ K it is fixed), with the
  extrapolation `1 + ln τ₁/(ln τ₁ + 1)` of the fit beyond τ = 1. Type 62 stops the polynomial three coefficients
  before the end and adds `c[m-2] ln(c[m-1] τ) e^{-c[m] τ}`. The third integer of the records and the level
  numbers `de` that `calt6062` computes from them are not used. The fit is not clamped at 0 (no record gives a
  negative rate), the temperature is at least ΔE/50k, the weights have 1e-16 added, and the levels are ordered by
  energy and must be in `1:nlev`. The `hot` fixtures reach 3×10⁹ K.
- `CollisionHelikeSat` (ucalc 73) evaluates the fit of `calt73` (Lotz-like terms with the exponential integrals
  `E₁`, `E₂`, `E₃`; the database only has the `r = 1` form, the `r = 2` one is ported but not checked). `ucalc`
  uses the first coefficient both as an energy in Rydberg (in the fit) and as a wavelength in Å (the temperature
  floor, the Boltzmann factor and the energy of the rates), so the Boltzmann factor is that of a transition of
  `12398.4/coeffs[1]` eV, which is not the energy difference of the levels; reproduced. The levels only order the
  two, and `Ry_K_sat` and `thermal_bohr` are the constants of the fit (`ucalc`: 1.578876e5 K and 5.46538e-11).
- `CollisionFe19` (ucalc 81): a constant Υ (negative values are 0), with the maxwellian coefficient and
  detailed balance like `ElectronImpact1`. `ucalc` orders the two levels by energy and accepts any level up to `nlev`
  (the continuum included); unlike type 56 it has no minimum energy difference.
- The dielectronic recombination rates (`DielecRecomb1` 7, `DielecRecombH` 22, `TotDielecRecomb` 39) return `frate` to the
  ground level of the next ion (`init = 1`, `final = 0`), without levels. The `T0` and `T1` of type 7 are in
  units of 10⁴ K, as `ucalc` uses them with the temperature in those units (the header said K). `DielecRecombH` (Storey)
  has five coefficients (`coeffs`; the header mentions a sixth) and is zero above 6×10⁴ K; `TotDielecRecomb` (Badnell)
  holds `C` and `T` term arrays (`T` in K, divided by 10⁴ here). Type 8 (`DielecRecomb2`, Arnaud and Raymond: four terms,
  energies in eV) has no records in `atdb.fits`: it is checked against `ucalc` with synthetic records
  (`test/reference/ucalc/synthetic.jl`, `type08.in`). Types 4 (`AtomicLine`) and 5 (`TwoPhotonRad`) have no records
  either and are not ported: `ucalc` reads their reals at positions (type 4: the fifth as an atomic mass; type 5: the fifth and sixth)
  that the structs do not have, which cannot be settled without data.
- `RadiativeFeDecay` (ucalc 82) takes the levels of the Fe UTA ions from two record types: the regular
  ones (type 6) and the UTA levels (`AtomicLevelFe`, type 83: energy and weight only, 986 records, for
  ions 327-341), so `levels` includes both, the second as an `AtomicLevel` with zero quantum numbers.
  The rate is that of `AtomicLine2` with these differences: the decay is `A_rad · pesc` (no floor, `A_auto` is not
  used), the oscillator strength is the stored `gf`, the line energy comes from the wavelength, there is no
  covering fraction and the energy of the decays is not returned (`ucalc` has only the photoexcitation's).
  Type 83 has no rate of its own.
- `CollisionAPED` (ucalc 92) had its integers in the wrong places (`[lower, upper, kind, Z, ionN]`, not
  `[1, lower, upper, Z, ionN]`) and its temperature limits in the grid; the record is now `transition`,
  `kind`, `T_min`, `T_max`, `T_grid` and `Υ`. `calc_maxwell_rates` has some thirty kinds of fit; the
  database has only kind 113 (1,302 records) and 116 (12), interpolated electron collision strengths with
  13 and 16 points (kinds `100 + n`), which are all `rate` supports (other kinds throw an error). The
  interpolation is log-log, gives 0 outside the first `n − 1` points (the last point is never used) and
  the rates are 0 outside `T_min`..`T_max`. The first level of the record is the lower one, whatever the
  level energies, and the Boltzmann factor uses `kB_keV` (`ucalc` has 8.617385e-8) and `ups_coeff`
  (8.629e-6). With the default temperatures of the fixtures only a quarter of the cases fall inside the
  table (`T_min` is 1.15e5 K); `type92hot` uses hotter ones.
- `CollisionIonize` (ucalc 95) keeps the scaled temperature grid `x_grid` and the strengths `ρ` apart (the old
  constructor put `T0` in `ρ`, and the header omits the grid). `ucalc` always takes level 1 as the initial
  level, and its statistical weight, whatever the level of the record; records of rate type 5 end in the
  continuum and those of rate type 15 (a duplicate set) in level 1. Below the first grid point `ucalc` reads
  the left end of the interval from before the table (`T0`, about 10⁴, and the last `x_grid`), which makes the
  interpolation come out at about ρ(x₁); reproduced. The Saha prefactor is `saha_ci` (2.08e-22 for `ucalc`).
  The factor `ln 2` is the single-precision `0.693147` there, and `kT_eV` of `ucalc_constants()` is now the
  single-precision `0.861707` too, since the Boltzmann factors amplify its 3e-8 rounding.
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
| `ChargeExHp` (10) | `level`, `ion`, `a`, `b`, `c`, `d`, `T_min`, `T_max`, `E_k`, `extra` |

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

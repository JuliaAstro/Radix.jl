# Using Radix

Radix solves the ionization, heating and radiative transfer of a photoionized slab of gas, as XSTAR does (`xstar`, 1-dimensional, one pass). It reads the atomic data
(`atdb.fits`) and the Compton table (`coheat.dat`) of an XSTAR installation (`$HEADAS/refdata/atdb.fits`, or `ftools/xstar/data/` of the source), and does not need XSTAR to run.

## A slab

```julia
using Radix

model = Model("atdb.fits", "coheat.dat")                  # the elements and processes (abundances = :xdef, the default table of XSTAR)

run = slab(model;
    density = 1e4,                # cm⁻³ (or `pressure = 1e-7` dyn cm⁻² for constant pressure)
    column = 1e21,                # cm⁻²: the zones go on until this hydrogen column
    logξ = 1.0,                   # log of L/(n r²) (of L/(4π c P r²) for the pressure)
    luminosity = 1e46,            # erg s⁻¹ between 13.6 eV and 13.6 keV
    α = -1.0,                     # the source is E^α (or `spectrum = (E, L) -> blackbody(E, 0.5, L)`, `thermal`, `tabulated`)
    equilibrium = true,           # iterate the temperature of each zone to the thermal equilibrium (`false` keeps `T`)
    T = 10.0, xee = 1.0,          # the first guess, in 10⁴ K, and of the electron fraction
    emult = 1.0, taumax = 5.0, steps = 3,    # the thickness of the zones (XSTAR's `emult`, `taumax` and `nsteps`)
    vturb = 0.0)                  # km s⁻¹: the widths of the lines
```

`Model(atdb, coheat; abundances, multiplier, exclude, cfrac)` takes an abundance table (`composition(:angr)`, ... or a vector for the elements 1 to 30),
`multiplier = Dict(26 => 2.0)` for factors, `exclude` for the elements that are left out (Li, Be and B by default, as the atomic data of XSTAR has no data for them) and the covering fraction.
Under the hood `slab` is `powerlaw` (the source), `source_distance`, `Base.iterate` of a mixture (the zones and the transfer); they can be called on their own, and it takes the zones (`(r, Δr)`) as well.

## What comes back

A named tuple:

- `zones`: one named tuple per zone, with `r` and `Δr` (cm), `T` (10⁴ K), `xee` and `nₑ`, `ntot` (the hydrogen density), `heating`, `cooling` and `imbalance` (the relative difference), `fractions[k][i]` (the fraction of ion `i` of the element `k` of
  `model.mixture.Z`), `populations` (of the levels), `radiation` (the field of the zone), `opacity` (`total`, `continuum`, `emissivity`, `bremsstrahlung`) and `escape`. The first zone has no thickness: it only gives the opacity from which the others are chosen.
- `E` (eV) and `incident` (erg s⁻¹ erg⁻¹), the source on XSTAR's energy grid, `transmitted` (`incident × exp(-dpthcont)`, the column of the output of XSTAR) and `r`;
- `spectrum`: the spectrum after the slab (`L`), the diffuse emission of the zones (`inward`, `outward`, and with the continuum opacity `inward_continuum`, `outward_continuum`, the columns `emit_inward` and `emit_outward` of XSTAR, in erg s⁻¹ erg⁻¹ like the incident spectrum);
- `luminosities` (`inward`, `outward`): the luminosity that each line and recombination edge of each element adds (10³⁸ erg s⁻¹; the records are `model.mixture.elements[k].rates`);
- `depths` (the optical depths of the lines and edges), `dpthc` and `dpthcont` (the depth of the continuum with and without lines and free-free).

Temperatures are in units of 10⁴ K, energies in eV, as in XSTAR; the spectra are in erg s⁻¹ erg⁻¹ (XSTAR's files are in units of 10³⁸ erg s⁻¹).

## Results as arrays, and one file

XSTAR writes a dozen files (`xout_abund1`, `xout_cont1`, `xout_spect1`, `xout_lines1`, `xout_rrc1`, the `xo_detal` files, a log, ...). Radix returns everything in memory, as the named tuple above, so that it can be called repeatedly
as a model (a fit, for instance) without touching the disk. `tables` turns the run into `Results`, four tables of plain arrays, one for each of the essential outputs:

```julia
t = tables(model, run)
t.zones       # radius, depth, thickness, ion_parameter, x_e, n_h, n_e, pressure, temperature, heating, cooling, heat_error and one column per ion (h_i, h_ii, he_i, ...)
t.spectrum    # energy, incident, transmitted, transmitted_lines, emit_inward, emit_outward, diffuse_inward, diffuse_outward, lines_inward, lines_outward (erg/s/erg, on XSTAR's energy grid)
t.lines       # ion, lower, upper, wavelength, emit_inward, emit_outward (erg/s), depth_inward, depth_outward
t.edges       # ion, level, energy, emit_inward, emit_outward (erg/s), depth_inward, depth_outward
```

`write("model.fits", t; meta=(; DENSITY=1e4, COLUMN=1e21))` writes the four tables to one FITS file (an empty primary HDU with the `meta` keywords and the extensions `ZONES`, `SPECTRUM`, `LINES` and `EDGES`, with units); it
passes HEASoft's `ftverify`. Nothing is written unless it is asked for. The spectrum has the pieces of both spectral files of XSTAR: `xout_cont1` has `transmitted`, `emit_inward` and `emit_outward` (the continuum), and `xout_spect1` has `transmitted_lines` and the emission `diffuse + lines` (the lines binned with their Voigt profiles, `binned=false` skips them: they take 10 s). Not in the tables: the detailed level populations (`run.zones[k].populations`
keeps them in memory).

## Agreement with XSTAR

Radix has been compared with XSTAR 2.59j for the rates of every record type, the ionization balance, the heating and cooling, the transfer (escape, optical depths, continuum opacity, lines, two-photon continua, diffuse emission), the thickness of the zones, the sources and constant pressure
([the field notes](field-names.md) have the details and the numbers, `test/reference/` the reference output). For the thick slab at log ξ = 1 it reproduces the zones of XSTAR (depths to 10⁻⁴), the temperatures to 0.05%, the ion fractions to 10⁻³ or better and the diffuse emission to 5%.
What differs:

- XSTAR ends its temperature iteration early in some zones (its `dsec` stops when the false position repeats a temperature): Radix finds the equilibrium, which is 0.4% to 2% from XSTAR's temperature in those zones.
- The atomic database of the XSTAR package (`$HEADAS/refdata/atdb.fits`) is not the one of its source tree: hydrogen has other records, and the temperatures of the tree are 3-6% lower than the package's. Use the one that you want to compare with.
- With `vturbi > 0` XSTAR 2.59j smooths its continuum so that it is not absorbed below 20 keV (`gsmooth2`): Radix does not, and the comparisons with it use `vturbi = 0`.
- To reproduce the rounded constants of XSTAR (`0.861707 eV` for k × 10⁴ K, 12.56 for 4π, ...) use `constants!(ucalc())` first; the default is CODATA 2022.
- The covering fraction `cfrac` cuts the photoexcitation of the lines, the escape probabilities and Thomson scattering, as in XSTAR (it agrees to 0.1% in temperature for `cfrac = 0.5`).
- The blackbody and thermal bremsstrahlung sources agree with XSTAR with the temperature iterated: the first zone of a blackbody (kT = 431 eV) to 0.03% in T and 0.3% in the ion fractions, and of bremsstrahlung (kT = 86 eV) to 0.4%. Through a slab of 10¹⁹ cm⁻² the first four zones of the blackbody (T from 11.8 to 7.6×10⁴ K) and all of the bremsstrahlung agree to 0.4% in T and 1% in the fractions. The blackbody slab then reaches a thermal front, where T falls to 1.8×10⁴ K in one zone: T agrees to 0.3% in the first zone of it and 1.7% at the end, and the trace ions of C to Fe differ by factors. XSTAR is not converged there (its T at the same radius changes by 1% with the size of the steps, it needs 20-30 iterations per zone, and its final row has a heating imbalance of 0.6%), and Radix changes by 0.8% in the same way.
- At a density of 10⁹ cm⁻³ the temperature of a zone at log ξ = 2 is 1.2% below XSTAR's (19.74 against 19.98, in the zone where the heating and cooling of XSTAR balance to 5×10⁻⁵); at 10⁴ cm⁻³ the same case agrees to 0.02%. It is the aluminium, and in it the helium-like Al XII: the heating minus the cooling of every other element agrees with XSTAR's table to a few 10⁻⁹. The rates of all the records of aluminium agree with XSTAR's `ucalc` at that cell (`test/reference/ucalc/generate.jl` with `CONDITION=` samples them), the photoexcitation rates agree to 0.02%, and the matrix of Radix has the same entries as the one that `msolvelucy` prints at `lprint=5` (both count the identical pairs of `CollisionAPED` records of Al XII that the package database has, and the 36 stale rates that `ucalc` leaves between levels of equal energy change nothing). What differs are the populations: XSTAR's 2³P levels are 8% lower than Radix's (4.81×10⁻¹¹ against 5.20×10⁻¹¹ for J = 0), so the three lines of the 2³S–2³P triplet (900, 943 and 954 Å), which are 80% of the line cooling of aluminium, are 8% weaker. XSTAR's populations are not a steady state of its own matrix: one sweep of the relaxation that `msolvelucy` runs after the solution of its condensed (one level for each ion) matrix, from XSTAR's final populations, changes the 2³P levels by +10%. That relaxation stops when the populations above 10⁻⁶ change by less than 1% (`crit2 = 10⁻²`), so the small populations that depend on a slowly converging level are left a sweep behind. Radix solves the whole matrix, and its populations satisfy it to 10⁻¹³.
- The emission of the thick slab agrees with XSTAR's to 3% in all the bins of its spectrum but one with the package database (the one is at 1.145 eV). That bin and the 5876 Å bin hold the triplets of He I, which XSTAR gives 1.8 to 2.2 times weaker than Radix: its populations of the He I (all below the 10⁻⁶ that its relaxation looks at) are not a steady state of its own matrix, which Radix reproduces to 10⁻³, and the 2³S gets 1.5 times what it loses there (`test/reference/xstar_hhe`, the slab with hydrogen and helium only, where the table of XSTAR lists every line).
- `passes=3` (XSTAR's `npass`) repeats the march with the optical depths of the lines and edges between each zone and the outer edge, as the escape probabilities of the outer side. It changes the ion fractions by up to 2.5% on the thick slab and follows XSTAR's change to 0.15%, the temperatures of the zones (with the temperature iterated) by 0.9%, and the lines much more (the He II 3s-2p line by a factor 1.6, the Fe X 3p-3d line by 15), which Radix follows: of the 566 lines of XSTAR's table that both have in the two runs, 545 change as in XSTAR to 2% (two differ by more than 5%, the forbidden C V line at 41 Å, +10% in XSTAR and 0% in Radix, and a line of Ar IX, +16% and +5%, because XSTAR's populations of their upper levels are below the 10⁻⁶ that its relaxation looks at and are not a steady state of its own matrix: `test/reference/xstar_hcar`). The continuum emission agrees to a median of 1.001 and 2.6% at worst, the transmitted spectrum to 3×10⁻⁴, the luminosities of the lines to a median of 1.001 (98% within 5%) and the total of the recombination edges to 0.2% (`test/reference/xstar_thick_slab_npass3`). `passes` must be odd.
- **Recombination continua narrower than a bin** (`exact_binning` of the constants): the recombination continuum above a K edge falls in about kT (9 eV at 10⁵ K) and is one or two bins wide at 999 bins (21 eV at 1.4 keV). XSTAR puts it in the continuum arrays as the value at the grid points, starting at the bin nearest to the threshold, so the integral of the emission of an edge depends on where its threshold falls in a bin: for the ground edge of Ne X it is 0.26, 1.65, 1.22 and 0.89 of the converged one at 999, 1500, 2000 and 3000 bins and 0.90 at 10000 (O VII 1.19, 1.29, 1.33, 0.89 and 0.96, O VIII 0.64 to 1.20), and 1.00 only at 30000 or more. The soft edges (He II at 54 eV, ...) are wide against a bin and agree to 1% at any number. The total emission of the thick slab does not show it (it averages over the edges), the ranges of 1-20 keV do: ±35% in the integral between 999 and 30000 bins. Radix does that with the constants of `ucalc()` (what the comparisons with XSTAR use) and by default takes the average of the opacity and of the emission of the edge over the cell of each bin (from the middle of the interval to the bin before to the middle of that to the one after, the weights of the trapezoid rule that integrates the arrays), which conserves the absorption and the energy of an edge whatever the position of the threshold: the integral of each of the 42 edges tested is the same at 999 and at 100000 bins to 4×10⁻⁷ and agrees with the converged value of XSTAR's rule to 0.2%. The balance (rates, heating and cooling) is still integrated with XSTAR's rule on its 999 bins.
- **The edges and lines per bin** (`rank_per_bin` of the constants, `nrank` of `iterate`, `slab` and `tables`): XSTAR puts only the ten strongest edges and lines of each energy bin in the continuum arrays and the binned lines, so what it calculates depends on the number of bins (the total continuum emission of the thick slab changes by 13% from 999 to 30000 bins, and by a factor 2 in 100-300 eV). Radix keeps (nearly) all of them by default (`rank_per_bin = 1000`), which with the exact binning of the edges below (`exact_binning`) makes the total continuum emission of the thick slab the same to 0.7% from 999 to 30000 bins (5.89, 5.85, 5.88 and 5.93 ×10¹⁷, against 5.83, 5.75, 5.83 and 5.93 with XSTAR's rule) and the energy of each edge the same to 10⁻⁶ (the integrals over a range of energy, which a sharp edge crosses at the boundary or not, still wander by a few percent at 999 bins: they depend on where the boundary falls among the grid points); `Radix.constants!(Radix.ucalc())` (the constants of XSTAR's source, with its 10) is what the comparisons with XSTAR of the tests use.
- **The number of energy bins** (`ncn2`: 10000 by default in `slab` and in `xstar_energy_grid`, 999 in the references of XSTAR, which the tests pass as `ncn2=999`; XSTAR's own default is 9999): XSTAR solves the ionization balance, the heating and cooling, and the edges on the field mapped to 999 bins whatever `ncn2` is (`bremsmap`), and only the continuum arrays (opacity and recombination emission, the transmitted spectrum, the binned lines) are on the `ncn2` bins; Radix does the same. The thin slab then changes by 0.03% in temperature from 999 to 10000 bins in XSTAR (18.8925, 18.8982) and by 0.03% in Radix (18.8958, 18.9015), and at 10000 bins Radix agrees with XSTAR to 0.02% in temperature and 0.3% in the ion fractions; the thick slab (with XSTAR's zones) to a median of 1.001 in the continuum emission (5% at worst, 10 of 5700 bins beyond 3%), 3×10⁻⁴ in the transmitted spectrum, and 0.1-0.8% in the fractions (`test/reference/xstar_10000bins`). The thickness of the zones of the thick slab with `step` is not the same after the eighth zone at 10000 bins: the rule takes the strongest opacity of the bins with a depth below `taumax = 5`, and the bin of the O VII line at 575 eV is at 4.988 in XSTAR and 5.025 in Radix (0.75% in the line), so XSTAR's zone is 4.2×10¹⁵ cm thick and Radix's 7.7×10¹⁵ (the 15 rows of XSTAR's table against 15 zones of Radix at other depths: the sums of the emission do not depend on the partition). With 999 bins the lines are averaged over wider bins and no depth is near the threshold. With exact binning the total continuum emission of the thick slab, and the energy of every edge, are independent of the number of bins (0.7% in the total from 999 to 30000 bins), and the absorption of the continuum edges is within 2% of the converged value from 3000 bins (1-14% at 999). The absorption by the lines is not either with XSTAR's mean opacity of the bins: the lines are narrower than a bin at any number (the thermal width is 6×10⁻⁵ of the energy, the bin 1.6×10⁻² at 999 bins and 1.6×10⁻³ at 10000), the transmission of a bin `exp(-⟨τ⟩)` is not the average `⟨exp(-τ)⟩` of a saturated line, and the equivalent width of the absorption of the thick slab with its lines is 1.5-3.4 times that of 30000 bins at 999 bins and 1.1-1.6 at 10000. The transmitted spectrum with the lines (`transmission` of `slab`, `transmitted_lines` of `tables`) takes by default (`resolved_lines` of the constants, `false` in `ucalc()`) the profile of each strong line inside the cells of the bins (`resolved_transmission`: the equivalent width of the saturated lines, the mean depth of the others), and then the equivalent width is within 7% of the converged one at 999 bins, 2% at 3000 and 1% from 10000 bins. The balance does not depend on it (it is on 999 bins as in XSTAR), so the cost does not either: 53 s at 999 bins and 62 s at 10000 bins for the thick slab.

## Speed

Radix runs the elements of the gas in parallel, so start Julia with several threads (`julia -t auto`): the default is one. Wall times on an Apple M-series laptop with the database of the XSTAR source tree (the first call also compiles):

| model | XSTAR 2.59j | Radix, 1 thread | Radix, 8 threads |
|---|---|---|---|
| thick slab at log ξ = 1, T and x_e kept, 7 zones, with the emission of all lines and edges | 148 s | 51 s | 28 s |
| slab of 10¹⁹ cm⁻² at log ξ = 2, T iterated | 411 s (10 zones) | 137 s (7 zones) | 78 s (7 zones) |

The zones are not the same in the second row (the thickness rule of Radix gives fewer), so compare the time of a zone: 41 s for XSTAR, 20 s and 11 s for Radix. `passes=3` takes about twice as long as `passes=1`, XSTAR's `npass=3` four times.
Most of the time is the integration of the photoionization cross sections over the radiation (once for each balance of the elements, the heating, the opacities and the emissivities of each zone).

## Running the tests

`Pkg.test("Radix")` runs the tests of the rates and the toy models. The comparisons with XSTAR (the level balance, heating, transfer, spectra and the whole database) need its atomic database, which is too large for the repository (875 MB): give the path of `atdb.fits` in `RADIX_ATDB` (the one of the XSTAR 2.59j package, `$HEADAS/refdata/atdb.fits`, or of its source tree: both pass) and use several threads:

```
RADIX_ATDB=$HEADAS/refdata/atdb.fits julia -t 4 --project -e 'using Pkg; Pkg.test()'
```

It takes about 7 minutes with 8 threads. `RADIX_QUICK=1` skips the five tests that take minutes (the thermal equilibrium of zones, the thick slab with its emission lines, the three passes with the temperature iterated and the model loaded from the files): the rest take about 3 minutes. The continuous integration runs that quick suite: it downloads the database from HEASARC once (`https://heasarc.gsfc.nasa.gov/FTP/software/lheasoft/xstar/atdb.fits`) and keeps it in the cache of the workflow, and the whole suite is run locally.

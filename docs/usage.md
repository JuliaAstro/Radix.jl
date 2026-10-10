# Using Radix

Radix solves the ionization, heating and radiative transfer of a photoionized slab of gas, as XSTAR does (`xstar`, 1-dimensional, one pass). It reads the atomic data
(`atdb.fits`) and the Compton table (`coheat.dat`) of an XSTAR installation (`$HEADAS/refdata/atdb.fits`, or `ftools/xstar/data/` of the source), and does not need XSTAR to run.

## A slab

```julia
using Radix

model = Model("atdb.fits", "coheat.dat")                  # the elements and processes (abundances = :xdef, the default table of XSTAR)

run = slab_model(model;
    density = 1e4,                # cm⁻³ (or `pressure = 1e-7` dyn cm⁻² for constant pressure)
    column = 1e21,                # cm⁻²: the zones go on until this hydrogen column
    logξ = 1.0,                   # log of L/(n r²) (of L/(4π c P r²) for the pressure)
    luminosity = 1e46,            # erg s⁻¹ between 13.6 eV and 13.6 keV
    α = -1.0,                     # the source is E^α (or `spectrum = (E, L) -> blackbody(E, 0.5, L)`, `thermal_bremsstrahlung`, `tabulated`)
    equilibrium = true,           # iterate the temperature of each zone to the thermal equilibrium (`false` keeps `T`)
    T = 10.0, xee = 1.0,          # the first guess, in 10⁴ K, and of the electron fraction
    emult = 1.0, taumax = 5.0, steps = 3,    # the thickness of the zones (XSTAR's `emult`, `taumax` and `nsteps`)
    vturb = 0.0)                  # km s⁻¹: the widths of the lines
```

`Model(atdb, coheat; abundances, multiplier, exclude, cfrac)` takes an abundance table (`abundance_table(:angr)`, ... or a vector for the elements 1 to 30),
`multiplier = Dict(26 => 2.0)` for factors, `exclude` for the elements that are left out (Li, Be and B by default, as the atomic data of XSTAR has no data for them) and the covering fraction.
Under the hood `slab_model` is `power_law` (the source), `source_distance`, `march_slab` (the zones) and `march_zones` (the transfer); they can be called on their own, and `march_zones` takes the zones (`(r, Δr)`) as well.

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
(`docs/field-names.md` has the details and the numbers, `test/reference/` the reference output). For the thick slab at log ξ = 1 it reproduces the zones of XSTAR (depths to 10⁻⁴), the temperatures to 0.05%, the ion fractions to 10⁻³ or better and the diffuse emission to 5%.
What differs:

- XSTAR ends its temperature iteration early in some zones (its `dsec` stops when the false position repeats a temperature): Radix finds the equilibrium, which is 0.4% to 2% from XSTAR's temperature in those zones.
- The atomic database of the XSTAR package (`$HEADAS/refdata/atdb.fits`) is not the one of its source tree: hydrogen has other records, and the temperatures of the tree are 3-6% lower than the package's. Use the one that you want to compare with.
- With `vturbi > 0` XSTAR 2.59j smooths its continuum so that it is not absorbed below 20 keV (`gsmooth2`): Radix does not, and the comparisons with it use `vturbi = 0`.
- To reproduce the rounded constants of XSTAR (`0.861707 eV` for k × 10⁴ K, 12.56 for 4π, ...) use `set_constants!(ucalc_constants())` first; the default is CODATA 2022.
- `passes=3` (XSTAR's `npass`) repeats the march with the optical depths of the lines and edges between each zone and the outer edge, as the escape probabilities of the outer side. It changes the ion fractions by up to 2.5% on the thick slab and follows XSTAR's change to 0.15%. `passes` must be odd.

## Speed

Radix runs the elements of the gas in parallel, so start Julia with several threads (`julia -t auto`): the default is one. Wall times on an Apple M-series laptop with the database of the XSTAR source tree (the first call also compiles):

| model | XSTAR 2.59j | Radix, 1 thread | Radix, 8 threads |
|---|---|---|---|
| thick slab at log ξ = 1, T and x_e kept, 7 zones, with the emission of all lines and edges | 148 s | 51 s | 28 s |
| slab of 10¹⁹ cm⁻² at log ξ = 2, T iterated | 411 s (10 zones) | 137 s (7 zones) | 78 s (7 zones) |

The zones are not the same in the second row (the thickness rule of Radix gives fewer), so compare the time of a zone: 41 s for XSTAR, 20 s and 11 s for Radix. `passes=3` takes about twice as long as `passes=1`, XSTAR's `npass=3` four times.
Most of the time is the integration of the photoionization cross sections over the radiation (once for each balance of the elements, the heating, the opacities and the emissivities of each zone).

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

## Agreement with XSTAR

Radix has been compared with XSTAR 2.59j for the rates of every record type, the ionization balance, the heating and cooling, the transfer (escape, optical depths, continuum opacity, lines, two-photon continua, diffuse emission), the thickness of the zones, the sources and constant pressure
(`docs/field-names.md` has the details and the numbers, `test/reference/` the reference output). For the thick slab at log ξ = 1 it reproduces the zones of XSTAR (depths to 10⁻⁴), the temperatures to 0.05%, the ion fractions to 10⁻³ or better and the diffuse emission to 5%.
What differs:

- XSTAR ends its temperature iteration early in some zones (its `dsec` stops when the false position repeats a temperature): Radix finds the equilibrium, which is 0.4% to 2% from XSTAR's temperature in those zones.
- The atomic database of the XSTAR package (`$HEADAS/refdata/atdb.fits`) is not the one of its source tree: hydrogen has other records, and the temperatures of the tree are 3-6% lower than the package's. Use the one that you want to compare with.
- With `vturbi > 0` XSTAR 2.59j smooths its continuum so that it is not absorbed below 20 keV (`gsmooth2`): Radix does not, and the comparisons with it use `vturbi = 0`.
- To reproduce the rounded constants of XSTAR (`0.861707 eV` for k × 10⁴ K, 12.56 for 4π, ...) use `set_constants!(ucalc_constants())` first; the default is CODATA 2022.
- Not done: more than one pass (`npass`), the density as a power of the radius (`radexp`), the spectra of the output files of XSTAR (the tables of lines and edges are in `luminosities` and `depths`), and the `delea` of the records of type 41 for the widths of the lines.

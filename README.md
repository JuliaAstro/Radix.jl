# Radix.jl

A spectral simulator for photo-ionized plasmas at X-ray frequencies

This package is based on Tim Kallman's XSTAR program that is written in FORTRAN. We believe that a Julia version will be more modular, extendable, and performant.

## Quick start

```julia
using Radix
model = Model("atdb.fits", "coheat.dat")        # the atomic data and the Compton table of an XSTAR installation
run = slab_model(model; density=1e4, column=1e21, logξ=1.0, luminosity=1e46, equilibrium=true, T=10.0, emult=1.0, steps=3, vturb=0.0)
[(z.T, z.xee, z.Δr) for z in run.zones]         # the zones: temperature (10⁴ K), electron fraction and thickness (cm)
```

See [docs/usage.md](docs/usage.md) for the options, the results and how Radix compares with XSTAR.


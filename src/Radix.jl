module Radix

using LinearAlgebra: lu, SingularException
using SparseArrays

export Atom, AtomicLevel, AtomicLine
export load, levels, counts, Levels, nlevels, weight, Radiation, Opacity
export Constants, ucalc, constants, constants!
export Model, slab, tables, march, ConstantPressure, composition
export powerlaw, blackbody, thermal, tabulated


include("constants.jl")
include("physical.jl")
include("atom.jl")
include("cell.jl")
include("ion.jl")

include("transition.jl")
include("chianti.jl")
include("hydrogenic.jl")
include("radiation.jl")
include("spectrum.jl")
include("photoionization.jl")
include("supertable.jl")
include("abstractrate.jl")
include("levels.jl")
include("rates/rates.jl")
include("ratemap.jl")

include("atomicdb.jl")
include("elements.jl")
include("heating.jl")
include("lines.jl")
include("processes.jl")
include("mixture.jl")
include("transfer.jl")
include("continuum.jl")

include("abundances.jl")
include("model.jl")
include("results.jl")

end

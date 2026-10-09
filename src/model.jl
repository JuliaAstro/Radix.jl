# The model of a run: the atomic data, the abundances and the processes of the gas, from the files of XSTAR.

"""
    Model(atdb, coheat; abundances=:xdef, multiplier=Dict(), exclude=(3, 4, 5), cfrac=0.0)

The gas of a run: the elements of the atomic database file `atdb` (`atdb.fits`) with the abundance table `abundances` (a name of `abundance_table`, or a table of the elements 1 to 30), their `multiplier`s (`Dict(Z => factor)`, XSTAR's
`habund`, `heabund`, ...), without the elements `exclude`d (Li, Be and B, which the atomic data of XSTAR does not have), and the processes of the gas that do not depend on the levels (`standard_processes`) with the Compton table `coheat`
(`coheat.dat`) and the covering fraction `cfrac`. `Model(records, compton; ...)` takes the loaded database and the `Compton` table. `slab_model(model; ...)` runs it and `tables(model, run)` makes its results into arrays. The fields are `mixture`, `processes` and `ion_labels`.
"""
struct Model{M<:Mixture, P}
    mixture::M
    processes::P
    ion_labels::Dict{Int, String}       # the labels of the ions of the database (`h_i`, `he_ii`, ...), by ion index
end

function Model(records::AbstractVector, compton::Compton; abundances=:xdef, multiplier=Dict{Int, Float64}(), exclude=(3, 4, 5), cfrac=0.0)
    factors = merge(Dict(Int(Z) => 0.0 for Z in exclude), Dict{Int, Float64}(multiplier))
    table = abundances isa Symbol ? abundance_table(abundances) : abundances
    Model(Mixture(records, levels(records); abundances=table, multiplier=factors), standard_processes(compton; cfrac),
        Dict{Int, String}(Int(r.ion) => String(strip(r.label)) for r in records if r isa Ion))
end

Model(atdb::AbstractString, coheat::AbstractString; kw...) = Model(open(load, atdb), load(Compton, coheat); kw...)

slab_model(model::Model; kw...) = slab_model(model.mixture, model.processes; kw...)

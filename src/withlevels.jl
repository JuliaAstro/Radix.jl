# Rates that need atomic data besides their own coefficients (level energies and weights, the number of levels of
# the ion, the atomic mass of its element) carry a `LevelTable` in a field `levels`.

"""
    @with_levels struct Name{...} <: AbstractRate ... end

Defines the struct with an extra last field `levels::L` (a new type parameter `L`) that holds the `LevelTable` the
rate reads its level data from, and the constructor with the original fields that leaves it `nothing`. Fill it with
`attach_levels`; `load` does so for every record.
"""
macro with_levels(ex)
    (ex isa Expr && ex.head === :struct) || error("@with_levels expects a struct definition")
    head, body = ex.args[2], ex.args[3]
    decl, sup = head isa Expr && head.head === :<: ? (head.args[1], head.args[2]) : (head, nothing)
    name, params = decl isa Expr && decl.head === :curly ? (decl.args[1], decl.args[2:end]) : (decl, Any[])
    nfields = count(a -> a isa Symbol || (a isa Expr && a.head === :(::)), body.args)
    push!(body.args, :(levels::L))
    newdecl = Expr(:curly, name, params..., :L)
    ex.args[2] = sup === nothing ? newdecl : Expr(:<:, newdecl, sup)
    quote
        $(esc(ex))
        $(esc(name))(args::Vararg{Any, $nfields}) = $(esc(name))(args..., nothing)
    end
end

"""
    attach_levels(coef, levels::LevelTable)

A copy of the rate `coef` that holds `levels` (a rate without a `levels` field is returned as it is). `coef.levels`
is what the rates that need level data read.
"""
function attach_levels(coef::T, table::LevelTable) where {T<:AbstractRate}
    hasfield(T, :levels) || return coef
    vals = ntuple(i -> fieldname(T, i) === :levels ? table : getfield(coef, i), fieldcount(T))
    Base.typename(T).wrapper(vals...)
end
attach_levels(coefs::AbstractArray, table::LevelTable) = map(c -> c isa AbstractRate ? attach_levels(c, table) : c, coefs)

# the level table of a rate, with a clear error if none was attached
function levels_of(coef)
    coef.levels === nothing && throw(ArgumentError(
        "$(typeof(coef).name.name) needs level data: use `attach_levels(coef, level_table(records))` " *
        "(the records of `load` have it)"))
    coef.levels
end

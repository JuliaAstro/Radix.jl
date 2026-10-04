# The atomic data the rates need besides their own coefficients.

"""
    LevelTable

The atomic data the rates need besides their own coefficients: the levels of every ion (energies, weights and
quantum numbers, as `AtomicLevel`s keyed by `(ion, level)`), the number of levels of each ion (`nlevels`: the highest
index, the continuum level) and the atomic mass of the element of each ion (`atomic_mass`). It behaves as the
dictionary of levels. `load` builds one and passes it to the constructor of every rate that needs it (the field `levels`).
"""
struct LevelTable{K, V} <: AbstractDict{K, V}
    levels::Dict{K, V}
    nlev::Dict{Int, Int}        # ion => number of levels
    mass::Dict{Int, Float64}    # ion => atomic mass of its element (amu)
end

Base.length(t::LevelTable) = length(t.levels)
Base.iterate(t::LevelTable, args...) = iterate(t.levels, args...)
Base.getindex(t::LevelTable, key) = t.levels[key]
Base.get(t::LevelTable, key, default) = get(t.levels, key, default)
Base.haskey(t::LevelTable, key) = haskey(t.levels, key)
Base.keys(t::LevelTable) = keys(t.levels)
Base.values(t::LevelTable) = values(t.levels)

"""
    nlevels(levels, ion)

The number of levels of `ion` in a `LevelTable`: its highest level index, the continuum level.
"""
nlevels(levels::LevelTable, ion) = get(levels.nlev, Int(ion), 0)

"""
    atomic_mass(levels, ion)

The atomic mass (amu) of the element of `ion` in a `LevelTable`.
"""
function atomic_mass(levels::LevelTable, ion)
    haskey(levels.mass, Int(ion)) || throw(ArgumentError("no atomic mass for ion $ion in the level table"))
    levels.mass[Int(ion)]
end

"""
    level_counts(levels)

Dictionary `ion => number of levels` (the highest level index, the continuum
level) from a `level_table`.
"""
level_counts(levels::LevelTable) = copy(levels.nlev)

"""
    needs_levels(T)

Whether the rate type `T` has a field `levels`, the `LevelTable` it reads its level data from.
"""
needs_levels(T::Type) = :levels in fieldnames(Base.unwrap_unionall(T))

"""
    construct(T, rate, label, ivec, rvec, levels)

The record of rate type `T` (an entry of `ratemap`) from the integers `ivec` and reals `rvec` of the database,
with the `LevelTable` `levels` if the type needs it (`needs_levels`).
"""
construct(T::Type, rate, label, ivec, rvec, levels) =
    needs_levels(T) ? T(Int32(rate), label, ivec, rvec, levels) : T(Int32(rate), label, ivec, rvec)

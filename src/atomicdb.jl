####    Atomic Database    ####

using FITSFiles

function load(io::IO)

    # scale=false: FITSFiles.jl (before JuliaAstro/FITSFiles.jl#50 is released)
    # rounds unscaled Int32 columns through Float32, which corrupts the pointers
    # to the reals beyond 2^24. The database has no TSCAL/TZERO, so nothing is lost.
    atomdb = fits(io; scale=false)

    ptr  = reshape(atomdb[2].data[:pointers][1], 10, :)
    rdat = atomdb[3].data[:reals][1]
    idat = atomdb[4].data[:integers][1]
    sdat = atomdb[5].data[:char][1]

    N = length(ptr[1,:])
    nrec = something(findfirst(j -> ptr[2, j] == 0, 1:N), N + 1) - 1     # the records end at the first empty one

    # the record j, with the level data if its rate type needs it
    function record(j, levels)
        rate = ratemap[ptr[2,j]]
        label = String(sdat[ptr[10,j]:ptr[10,j]+ptr[7,j]-1])
        ivec  = idat[ptr[ 9,j]:ptr[ 9,j]+ptr[6,j]-1]
        rvec  = rdat[ptr[ 8,j]:ptr[ 8,j]+ptr[5,j]-1]
        construct(rate, ptr[3,j], label, ivec, rvec, levels)
    end
    implemented(j) = !(ratemap[ptr[2,j]] in (Missing, Nothing))

    rates = Array{Union{Atom, Ion, AbstractRate, Missing, Nothing}}(missing, N)
    for j in 1:nrec
        implemented(j) || println("Rate $(ptr[2,j]) not implemented.")
    end
    # the levels, the elements and the ions first: the other rates carry the level table built from them
    for j in 1:nrec
        implemented(j) && ptr[2,j] in level_data_types && (rates[j] = record(j, nothing))
    end
    table = level_table(rates)
    for j in 1:nrec
        implemented(j) && !(ptr[2,j] in level_data_types) && (rates[j] = record(j, table))
    end
    rates
end

const level_data_types = (6, 13, 14, 83)   # AtomicLevel, Atom, Ion and AtomicLevelFe records

"""
    level_table(records; masses=Dict())

The `LevelTable` of a loaded database (or any collection of records): the `AtomicLevel` records, the Fe UTA levels
(`AtomicLevelFe`, which hold only the energy and the weight, added as `AtomicLevel`s with zero quantum numbers unless
the level is already there), and the atomic masses of the ions from the `Atom` and `Ion` records. `masses`
(`ion => mass`) adds to or replaces those.
"""
function level_table(records; masses=Dict{Int, Float64}())
    levels = Dict((r.ion, r.level) => r for r in records if r isa AtomicLevel)
    for r in records
        r isa AtomicLevelFe || continue
        haskey(levels, (r.ion, r.level)) && continue
        z, zr = zero(r.level), zero(r.E)
        levels[(r.ion, r.level)] = AtomicLevel(r.rtype, r.label, z, z, z, z, r.level, r.ion, r.E, r.g, zr, zr)
    end
    nlev = Dict{Int, Int}()
    for (ion, level) in keys(levels)
        nlev[Int(ion)] = max(get(nlev, Int(ion), 0), Int(level))
    end
    element_mass = Dict(Int(r.Z) => Float64(r.mass) for r in records if r isa Atom)
    mass = Dict{Int, Float64}()
    for r in records
        r isa Ion && haskey(element_mass, Int(r.Z)) && (mass[Int(r.ion)] = element_mass[Int(r.Z)])
    end
    merge!(mass, Dict(Int(k) => Float64(v) for (k, v) in masses))
    LevelTable(levels, nlev, mass)
end


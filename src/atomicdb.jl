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

    rates = Array{Union{Atom, Ion, AbstractRate, Missing, Nothing}}(missing, N)
    for j=1:N
        if ptr[2,j] == 0
            break
        end
        rate = ratemap[ptr[2,j]]
        if rate in (Missing, Nothing)
            println("Rate $(ptr[2,j]) not implemented.")
            continue
        end
        type  = ptr[3,j]
        label = String(sdat[ptr[10,j]:ptr[10,j]+ptr[7,j]-1])
        ivec  = idat[ptr[ 9,j]:ptr[ 9,j]+ptr[6,j]-1]
        rvec  = rdat[ptr[ 8,j]:ptr[ 8,j]+ptr[5,j]-1]
        # if ptr[2,j] in (2, 9, 10, 13, 14)
        #     println("$(ptr[2,j]) => $rate: $ivec, $rvec, $label")
        # end
        rates[j] = rate(type, label, ivec, rvec)
    end
    rates
end

"""
    level_table(records)

Dictionary `(ion, level) => AtomicLevel` over the `AtomicLevel` records of a
loaded database, for the rates that need level energies and weights. The Fe UTA levels
(`AtomicLevelFe`, which hold only the energy and the weight) are added as `AtomicLevel`s
with zero quantum numbers, unless the level is already there.
"""
function level_table(records)
    table = Dict((r.ion, r.level) => r for r in records if r isa AtomicLevel)
    for r in records
        r isa AtomicLevelFe || continue
        haskey(table, (r.ion, r.level)) && continue
        z, zr = zero(r.level), zero(r.E)
        table[(r.ion, r.level)] = AtomicLevel(r.rtype, r.label, z, z, z, z, r.level, r.ion, r.E, r.g, zr, zr)
    end
    table
end

"""
    level_counts(levels)

Dictionary `ion => number of levels` (the highest level index, the continuum
level) from a `level_table`.
"""
function level_counts(levels)
    n = Dict{eltype(first(keys(levels))), Int}()
    for (ion, level) in keys(levels)
        n[ion] = max(get(n, ion, 0), Int(level))
    end
    n
end

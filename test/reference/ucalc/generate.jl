# Sample real records of one XSTAR data type from atdb.fits and write the input
# for the ucalc driver (drvu): one case per record and (t, nH, xe) condition,
# with the level data of the record's ion.
#
#   julia generate.jl /path/to/atdb.fits OUTDIR TYPE [NRECORDS]
#
# then:  drvu < OUTDIR/typeNN.in > OUTDIR/typeNN.out
using FITSFiles, Random

const CONDITIONS = (            # t [1e4 K], xpx [cm⁻³], xee, xh0/xpx, xh1/xpx
    (0.3, 1e4, 0.5, 0.5, 0.5),
    (1.0, 1e8, 0.9, 0.1, 0.9),
    (10.0, 1e10, 1.0, 1e-3, 1.0),
    (100.0, 1e12, 1.0, 1e-6, 1.0),
)

function main(path, outdir, type, nrec = 25)
    a = open(io -> fits(io; scale = false), path)
    ptr  = reshape(a[2].data[:pointers][1], 10, :)
    rdat = a[3].data[:reals][1]
    idat = a[4].data[:integers][1]
    N = findfirst(j -> ptr[2, j] == 0, 1:size(ptr, 2)) - 1
    ints(j)  = Int.(idat[ptr[9, j]:ptr[9, j] + ptr[6, j] - 1])
    reals(j) = Float64.(rdat[ptr[8, j]:ptr[8, j] + ptr[5, j] - 1])

    # levels per ion: ion => [(level, E, g, Einf, n, 2S+1, L)]
    levels = Dict{Int, Vector{NTuple{7, Float64}}}()
    for j in 1:N
        ptr[2, j] == 6 || continue
        i, r = ints(j), reals(j)                       # n, 2S+1, L, Z, level, ion ; E, g, ν, E∞
        push!(get!(levels, i[6], NTuple{7, Float64}[]),
              (i[5], r[1], r[2], r[4], i[1], i[2], i[3]))
    end
    foreach(v -> sort!(v, by = first), values(levels))

    # atomic masses: element records are [n_ions, Z], reals [abundance, mass]
    mass = Dict{Int, Float64}()
    for j in 1:N
        ptr[2, j] == 13 && (mass[ints(j)[2]] = reals(j)[2])
    end

    cand = [j for j in 1:N if ptr[2, j] == type]
    sel = cand[randperm(MersenneTwister(type), length(cand))[1:min(nrec, length(cand))]]
    mkpath(outdir)
    open(joinpath(outdir, "type$(lpad(type, 2, '0')).in"), "w") do io
        for j in sel
            iv, rv = ints(j), reals(j)
            lv = get(levels, iv[end], NTuple{7, Float64}[])
            for (t, xpx, xee, f0, f1) in CONDITIONS
                println(io, type, " ", ptr[3, j], " ", length(rv), " ", length(iv))
                println(io, join(iv, " "))
                println(io, join(rv, " "))
                amass = type == 50 ? get(mass, iv[3], 1.0) : 1.0       # type 50 ints: i, k, Z, ionN
                println(io, join((t, xpx, xee, xpx*f0, xpx*f1, 1.0, -1.0, 0.0, 1e-3, 1e-3, 0.5, 0.5, amass), " "))
                # nlev, then only the levels the record refers to
                need = sort(unique(filter(i -> 1 <= i <= length(lv), iv[1:end-1])))
                println(io, length(lv), " ", length(need))
                for i in need
                    (_, E, g, Einf, n, s2, L) = lv[i]
                    println(io, join((i, E, g, Einf, Int(n), Int(s2), Int(L)), " "))
                end
            end
        end
    end
    println("type $type: ", length(sel), " records x ", length(CONDITIONS), " conditions")
end

main(ARGS[1], ARGS[2], parse(Int, ARGS[3]), length(ARGS) > 3 ? parse(Int, ARGS[4]) : 25)

# Synthetic inputs for the ucalc driver (drvu) for rate types that have no record in atdb.fits,
# so that they can be checked against ucalc too (the formats are those of generate.jl).
#
#   julia synthetic.jl OUTDIR
#
# writes OUTDIR/type08.in (DielecRecomb2); then  drvu < OUTDIR/type08.in > OUTDIR/type08.out
using Random

const CONDITIONS = (            # t [1e4 K], xpx [cm⁻³], xee, xh0/xpx, xh1/xpx
    (0.3, 1e4, 0.5, 0.5, 0.5),
    (1.0, 1e8, 0.9, 0.1, 0.9),
    (10.0, 1e10, 1.0, 1e-3, 1.0),
    (100.0, 1e12, 1.0, 1e-6, 1.0),
)

# a case: the record (type, rtype, ints, reals), at every condition, with two levels (1 and the continuum 2)
function write_case(io, type, rtype, ints, reals)
    for (t, xpx, xee, f0, f1) in CONDITIONS
        println(io, type, " ", rtype, " ", length(reals), " ", length(ints))
        println(io, join(ints, " "))
        println(io, join(Float64.(reals), " "))      # (the exact doubles of the Float32 values, as in the database)
        println(io, join((t, xpx, xee, xpx*f0, xpx*f1, 1.0, -1.0, 0.0, 1e-3, 1e-3, 0.5, 0.5, 1.0), " "))
        println(io, join((0.1, 0.0015, 1.0, 0.0, 1, 1.0, 0.0, 1), " "))
        println(io, "2 2")
        println(io, "1 0.0 2.0 13.6 1 2 0")
        println(io, "2 13.6 1.0 13.6 1 2 0")
    end
end

function main(outdir)
    mkpath(outdir)
    rng = MersenneTwister(8)
    # type 8 (Arnaud and Raymond): four coefficients c (cm³ s⁻¹ K^{3/2}) and four energies E (eV)
    open(joinpath(outdir, "type08.in"), "w") do io
        for _ in 1:15
            c = Float32.(10 .^ (rand(rng, 4) .* 4 .- 6))
            e = Float32.(rand(rng, 4) .* 50)
            write_case(io, 8, 8, Int32[12], vcat(c, e))
        end
    end
    println("type 8: 15 records x ", length(CONDITIONS), " conditions")
end

main(ARGS[1])

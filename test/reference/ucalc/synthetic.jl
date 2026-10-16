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
function write_case(io, type, rtype, ints, reals; levels=((1, 0.0, 2.0), (2, 13.6, 1.0)), radiation=(0.1, 0.0015, 1.0, 0.0, 1, 1.0, 0.0, 1))
    for (t, xpx, xee, f0, f1) in CONDITIONS
        println(io, type, " ", rtype, " ", length(reals), " ", length(ints))
        println(io, join(ints, " "))
        println(io, join(Float64.(reals), " "))      # (the exact doubles of the Float32 values, as in the database)
        println(io, join((t, xpx, xee, xpx*f0, xpx*f1, 1.0, -1.0, 0.0, 1e-3, 1e-3, 0.5, 0.5, 1.0), " "))
        println(io, join(radiation, " "))
        println(io, length(levels), " ", length(levels))
        for (i, e, g) in levels
            println(io, i, " ", e, " ", g, " 13.6 1 2 0")
        end
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
    # type 65 (Sampson and Zhang's hydrogenic collisional ionization): n, L, 2J, Z, level, ion and the effective charge
    open(joinpath(outdir, "type65.in"), "w") do io
        for n in (1, 2, 3, 5, 8, 10, 11, 12, 15, 25), zeff in Float32.((1.0, 2.5, 8.0))
            write_case(io, 65, 1, Int32[n, 0, 1, 8, 1, 13], Float32[zeff])
        end
    end
    println("type 65: 30 records x ", length(CONDITIONS), " conditions")
    # type 67 (Keenan et al.): the levels (lower, upper), Z and the ion; the coefficients of the quadratic in log T (the first is the wavelength in Å)
    open(joinpath(outdir, "type67.in"), "w") do io
        for (levels, c1, c2, c3) in ((Int32[1, 2], 91.0, -12.0, 1.3), (Int32[2, 1], 91.0, -12.0, 1.3), (Int32[1, 2], 304.0, 40.0, -3.0), (Int32[1, 2], -304.0, 40.0, -3.0),
                (Int32[1, 2], 620.5, -250.0, 22.0), (Int32[1, 2], 1e-30, 1.0, 1.0), (Int32[2, 1], 12.4, 0.5, -0.04), (Int32[1, 2], 2000.0, -150.0, 1.0))
            write_case(io, 67, 1, Int32[levels..., 8, 13], Float32[c1, c2, c3])
        end
    end
    println("type 67: 8 records x ", length(CONDITIONS), " conditions")
    # type 96 (autoionization of satellite levels, Safronova): ionN, kN, ionN-1, iN-1, ionN and Aa, the factor, the energy above the limit (eV)
    open(joinpath(outdir, "type96.in"), "w") do io
        for (k, i, a, w, e) in ((1, 1, 1e12, 1.0, 0.0), (2, 1, 2e13, 3.0, 5.0), (1, 2, 5e12, 0.5, 30.0), (0, 3, 1e11, 2.0, 100.0), (2, 2, 3e12, 1.5, 1.2), (1, 1, 1e12, 1.0, -4.0))
            write_case(io, 96, 1, Int32[13, k, 12, i, 13], Float32[a, w, e])
        end
    end
    println("type 96: 6 records x ", length(CONDITIONS), " conditions")
    # type 79 (fluorescence lines): the levels, Z, ion; wavelength (Å), oscillator strength, two reals that ucalc does not read, the atomic mass
    open(joinpath(outdir, "type79.in"), "w") do io
        for (levels, λ, f, mass) in ((Int32[1, 2], 1000.0, 0.5, 16.0), (Int32[2, 1], 1000.0, 0.5, 16.0), (Int32[1, 2], 1.5, 0.02, 56.0), (Int32[1, 2], 304.0, 1e-13, 4.0),
                (Int32[1, 2], 2e9, 0.3, 1.0), (Int32[1, 2], -912.0, 0.1, 12.0), (Int32[1, 2], 6563.0, 0.64, 1.0))
            write_case(io, 79, 1, Int32[levels..., 8, 13], Float32[λ, f, 0.0, 0.0, mass])
        end
    end
    println("type 79: 7 records x ", length(CONDITIONS), " conditions")
    # type 89 (Safronova lines): the levels, Z, ion; wavelength (Å), a real that ucalc does not read, A (s⁻¹); the element has the atomic mass of the case (1)
    open(joinpath(outdir, "type89.in"), "w") do io
        for (levels, λ, a) in ((Int32[1, 2], 1000.0, 1e8), (Int32[2, 1], 1000.0, 1e8), (Int32[1, 2], 1.5, 3e12), (Int32[1, 2], 304.0, 1e10),
                (Int32[1, 2], 2e9, 5e5), (Int32[1, 2], -912.0, 2e9), (Int32[1, 2], 0.0, 1e9), (Int32[1, 2], 6563.0, 4.4e7), (Int32[1, 3], 100.0, 1e9), (Int32[3, 1], 100.0, 1e9), (Int32[1, 4], 100.0, 1e9))
            write_case(io, 89, 1, Int32[levels..., 8, 13], Float32[λ, 0.0, a]; levels=((1, 0.0, 2.0), (2, 10.2, 8.0), (3, 12.1, 18.0), (4, 13.6, 1.0)))
        end
    end
    println("type 89: 11 records x ", length(CONDITIONS), " conditions")
    # type 102 (SPEX collision strengths): itype, nreso, the two levels, Z, ion; the energy (keV), the six parameters of the continuum, four not used, the positions and the amplitudes of the resonances
    open(joinpath(outdir, "type102.in"), "w") do io
        for itype in 1:8, (nreso, levels) in ((0, (1, 2)), (1, (2, 1)), (2, (1, 3)), (3, (3, 2)), (4, (1, 2)), (4, (2, 4)))
            a = Float32.(0.1 .+ 2 .* rand(rng, 6))
            e = Float32(10 .^ (rand(rng)*2 - 3.5))          # 0.3 eV to 300 eV
            nres = nreso in (1, 2, 3) ? nreso : 4
            reso = Float32.(vcat(0.5 .+ 3 .* rand(rng, nres), 0.01 .+ 0.2 .* rand(rng, nres)))
            write_case(io, 102, 1, Int32[itype, nreso, levels..., 8, 13], vcat(e, a, zeros(Float32, 4), reso);
                levels=((1, 0.0, 2.0), (2, 10.2, 8.0), (3, 12.1, 18.0), (4, 13.6, 1.0)))
        end
    end
    println("type 102: 48 records x ", length(CONDITIONS), " conditions")
    # type 64 (hydrogenic photoionization): n, l, the charge, Z, the level, the ion. ucalc takes the top of the grid of the integral from `nphint`, which an earlier
    # call of type 59 leaves in its static memory, so that the file starts with a case of type 59 (the first of data/type59.in) that the test leaves out
    prime = readlines(joinpath(@__DIR__, "data", "type59.in"))[1:11]
    open(joinpath(outdir, "type64.in"), "w") do io
        println(io, join(prime, "\n"))
        for (n, l, charge, level) in ((1, 0, 1, 1), (2, 0, 1, 2), (2, 1, 1, 2), (3, 2, 1, 3), (3, 0, 2, 3), (1, 0, 2, 1), (2, 1, 3, 2), (3, 1, 1, 2))
            write_case(io, 64, 1, Int32[n, l, charge, 8, level, 13], Float32[0.0]; levels=((1, 0.0, 2.0), (2, 10.2, 8.0), (3, 12.1, 18.0), (4, 13.6, 1.0)),
                radiation=(0.1, 0.0015, 1.0, 1e16, 1, 1.0, 0.0, 1))
        end
    end
    println("type 64: a priming case and 8 records x ", length(CONDITIONS), " conditions")
end

main(ARGS[1])

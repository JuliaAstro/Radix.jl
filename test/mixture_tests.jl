# The ionization balance of a mixture of elements at the electron density that they produce together.

using FITSFiles

# ionization 1 → 2 at 1 s⁻¹ and recombination at α nₑ: the ion fraction of the continuum is 1/(1 + α nₑ)
struct RecombiningToyRate <: Radix.AbstractRate
    rtype::Int8
    ion::Int32
    α::Float64
end
Radix.rate(c::RecombiningToyRate, cell; kw...) = (; init=1, final=2, frate=1.0, irate=c.α*cell.nₑ)

function toy_mixture_tests()
    @testset "Mixture" begin
        f32 = Float32
        lv(level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1, 2, 0, 1, level, 5], f32[E, g, 1, 13.6])
        atom = Radix.Atom("hydrogen", Int32(1), Int32(1), 1.0f0, 1.01f0)
        helium = Radix.Atom("helium", Int32(2), Int32(2), 0.1f0, 4.0f0)
        ion = Radix.Ion("h_i", Int32(1), Int32(1), Int32(5), 13.6f0)
        records = [lv(1, 0.0, 2), lv(2, 13.6, 1), atom, helium, ion, RecombiningToyRate(Int8(3), Int32(5), 2e-4)]
        levels = Radix.levels(records)
        ntot = 1e4

        @testset "the elements" begin
            @test Radix.Mixture(records, levels; multiplier=Dict(2 => 0.0)).Z == [1]          # a zero abundance leaves an element out
            @test Radix.Mixture(records, levels; multiplier=Dict(2 => 1e-25)).Z == [1]        # as does one below `minimum_abundance`
            m = Radix.Mixture(records[[1, 2, 3, 5, 6]], levels)
            @test m.Z == [1] && m.abundance == [1.0] && length(m.elements) == 1
            @test Radix.Mixture(records[[1, 2, 3, 5, 6]], levels; multiplier=Dict(1 => 0.5)).abundance == [0.5]
        end

        @testset "the electrons" begin
            # the neutral atom has no electron free, an ion q, the bare nucleus Z; the abundance weights the elements
            m = Radix.Mixture(levels, [1, 2], [1.0, 0.1], Radix.Elements{typeof(levels)}[])
            @test Radix.electrons(m, [[0.2, 0.8], [0.1, 0.3, 0.6]]) ≈ 0.8 + 0.1*(0.3*1 + 0.6*2)
            @test Radix.electrons(m, [[1.0, 0.0], [1.0, 0.0, 0.0]]) == 0
            @test Radix.electrons(m, [[0.0, 1.0], [0.0, 0.0, 1.0]]) ≈ 1 + 0.1*2
        end

        @testset "the electron fraction" begin
            m = Radix.Mixture(records[[1, 2, 3, 5, 6]], levels)
            # the fraction of ionized hydrogen is x = 1/(1 + α ntot xee) and xee = x, so α ntot xee² + xee - 1 = 0 (α ntot = 2): xee = 1/2
            fixed = Radix.ionization_balance(m, 1.0, ntot; iterate=false)
            @test fixed.xee == 1 && fixed.nₑ == ntot && fixed.iterations == 1
            @test fixed.fractions[1][2] ≈ 1/3 && fixed.electrons ≈ 1/3                   # at xee = 1: α nₑ = 2
            @test fixed.converged                                                        # (nothing to iterate)
            solved = Radix.ionization_balance(m, 1.0, ntot)
            @test solved.converged
            @test solved.xee ≈ 0.5 rtol=2e-4
            @test solved.nₑ ≈ 5000 rtol=2e-4
            @test abs(solved.xee - solved.electrons)/solved.xee < 1e-4
            @test solved.fractions[1][2] ≈ solved.xee rtol=1e-4
            @test sum(solved.fractions[1]) ≈ 1 && sum(solved.populations[1]) ≈ 1
            # from either side, and with a tighter tolerance
            @test Radix.ionization_balance(m, 1.0, ntot; xee=3.0).xee ≈ 0.5 rtol=2e-4
            @test Radix.ionization_balance(m, 1.0, ntot; xee=0.01).xee ≈ 0.5 rtol=2e-4
            @test Radix.ionization_balance(m, 1.0, ntot; tolerance=1e-10).xee ≈ 0.5 rtol=1e-9
        end

        @testset "the bracketing of the iteration" begin
            # electrons 1.3 whatever the fraction: elcter = x - 1.3, found from below and from above (XSTAR's steps of 1.2, false position)
            for start in (1.0, 2.0, 0.1)
                xee, elcter, n, converged = Radix.bracket_electron_fraction(x -> x - 1.3, start, start - 1.3; tolerance=1e-4)
                @test converged
                @test xee ≈ 1.3 rtol=1e-4
                @test n <= 20
            end
            # not converging: the number of evaluations is limited
            xee, elcter, n, converged = Radix.bracket_electron_fraction(x -> x + 1.0, 1.0, 2.0; tolerance=1e-4)
            @test !converged && n == Radix.electron_iterations
        end
    end
end

# all the elements of the database: the electrons per hydrogen nucleus and the ion fractions of the XSTAR reference run
# (test/reference/xstar_pow_xi2), which keeps the electron fraction at 1
function mixture_balance_tests(db)
    @testset "Mixture of the database" begin
        levels = Radix.levels(db)
        mixture = Radix.Mixture(db, levels; multiplier=Dict(3 => 0.0, 4 => 0.0, 5 => 0.0))     # (the reference run has no Li, Be or B)
        @test length(mixture.Z) == 27 && !any(in(mixture.Z), (3, 4, 5))
        @test mixture.abundance[1:3] == [1.0, Float64(Float32(0.1)), Float64(Float32(0.00037))]

        dir = joinpath(@__DIR__, "reference", "xstar_pow_xi2")
        spectrum = fits(joinpath(dir, "xout_cont1.fits"))[3].data
        abundances = fits(joinpath(dir, "xout_abund1.fits"))[2].data
        roman = ["i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x", "xi", "xii", "xiii", "xiv", "xv", "xvi", "xvii", "xviii",
                 "xix", "xx", "xxi", "xxii", "xxiii", "xxiv", "xxv", "xxvi", "xxvii", "xxviii", "xxix", "xxx"]
        symbols = ["h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na", "mg", "al", "si", "p", "s", "cl", "ar", "k", "ca", "sc",
                   "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn"]
        row = 3                                             # the zone at 1.5×10¹³ cm (its outer edge is tabulated)
        reference = [vcat([getproperty(abundances, Symbol(symbols[Z], "_", roman[k]))[row] for k in 1:Z],
                          max(0.0, 1 - sum(getproperty(abundances, Symbol(symbols[Z], "_", roman[k]))[row] for k in 1:Z))) for Z in mixture.Z]
        radiation = Radix.point_source(Float64.(spectrum.energy), Float64.(spectrum.incident)*1e38, 1.5e13)

        @testset "at the electron fraction of the reference run" begin
            fixed = Radix.ionization_balance(mixture, 100.0, 1e4; radiation, iterate=false)
            @test fixed.electrons ≈ Radix.electrons(mixture, reference) rtol=1e-5
            worst = 0.0
            for (k, f) in enumerate(reference), i in eachindex(f)
                f[i] > 1e-3 && (worst = max(worst, abs(fixed.fractions[k][i]/f[i] - 1)))
            end
            @test worst < 5e-2                              # (every ion above 10⁻³ of every element)
        end

        @testset "the electron fraction of the gas" begin
            solved = Radix.ionization_balance(mixture, 100.0, 1e4; radiation)
            @test solved.converged && abs(solved.xee - solved.electrons)/solved.xee < 1e-4
            @test solved.xee ≈ Radix.electrons(mixture, reference) rtol=1e-3
            @test solved.nₑ == 1e4*solved.xee && 0 < solved.nₕ < 1e-2
            @test all(f -> sum(f) ≈ 1, solved.fractions)
        end
    end
end

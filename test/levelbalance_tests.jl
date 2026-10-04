# The level balance of one ion: the rate matrix and its populations.

using FITSFiles

# a rate between two levels of a toy ion that returns what it is told: (init, final, forward, inverse)
struct ToyRate <: Radix.AbstractRate
    rtype::Int8
    ion::Int32
    r::NamedTuple{(:init, :final, :frate, :irate), Tuple{Int, Int, Float64, Float64}}
end
Radix.rate(c::ToyRate, cell; kw...) = c.r

function toy_balance_tests()
  @testset "Level balance" begin
    f32 = Float32
    lv(level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1, 2, 0, 1, level, 5], f32[E, g, 1, 13.6])
    levels = Radix.level_table([lv(1, 0.0, 2), lv(2, 10.2, 8), lv(3, 12.0, 8), lv(4, 13.6, 1)]; masses=Dict(5 => 1.0))
    cell = Radix.Cell(1.0, 0.0, 1e4, 1e4)
    toy(i, f, a, b; type=3) = ToyRate(Int8(type), Int32(5), (; init=i, final=f, frate=a, irate=b))

    @testset "the matrix" begin
        # the energies order the two levels: the forward rate of (3, 2) is from level 3 (higher) to 2
        A = Radix.rate_matrix([toy(1, 2, 3.0, 5.0)], levels, 5, cell)
        @test A[2, 1] == 3.0 && A[1, 2] == 5.0 && A[1, 1] == -3.0 && A[2, 2] == -5.0
        # the first rate of a pair listed upper-first is still ucalc's rate lower -> upper
        A = Radix.rate_matrix([toy(2, 1, 3.0, 5.0)], levels, 5, cell)
        @test A[2, 1] == 3.0 && A[1, 2] == 5.0
        # photoionization keeps its order
        A = Radix.rate_matrix([toy(1, 4, 3.0, 5.0; type=7)], levels, 5, cell)
        @test A[4, 1] == 3.0 && A[1, 4] == 5.0 && A[1, 1] == -3.0 && A[4, 4] == -5.0
        A = Radix.rate_matrix([toy(4, 1, 3.0, 5.0; type=7)], levels, 5, cell)
        @test A[1, 4] == 3.0 && A[4, 1] == 5.0
        # records without a level pair and a level outside the ion do not enter, nor do the totals
        A = Radix.rate_matrix([toy(0, 0, 1.0, 1.0), toy(1, 0, 1.0, 1.0), toy(1, 9, 1.0, 1.0)], levels, 5, cell)
        @test all(iszero, A)
        @test isempty(Radix.ion_rates([toy(1, 2, 1.0, 1.0; type=8), toy(1, 2, 1.0, 1.0; type=15)], 5))
        @test length(Radix.ion_rates([toy(1, 2, 1.0, 1.0), ToyRate(Int8(3), Int32(6), (; init=1, final=2, frate=1.0, irate=1.0))], 5)) == 1
        # columns sum to zero: the transitions only move population
        A = Radix.rate_matrix([toy(1, 2, 3.0, 5.0), toy(2, 3, 1.0, 2.0), toy(1, 4, 7.0, 11.0; type=7)], levels, 5, cell)
        @test all(abs.(sum(A, dims=1)) .< 1e-12)
        # the line types of ucalc return the photoexcitation first; Radix names the decay the forward rate
        line = Radix.AtomicLine2(Int32(4), "", Int32[2, 1, 1, 5], f32[1215.67, 0.4162, 6.265e8], levels)
        A = Radix.rate_matrix([line], levels, 5, cell)
        @test A[1, 2] == Float64(f32(6.265e8)) && A[2, 1] == 0                      # decay 2 -> 1, no radiation
    end

    @testset "stacked ions" begin
        # ion 5 has the levels 1-3 and its continuum (4), ion 6 the level 1 and its continuum (2): the unknowns are
        # 5.1-5.3 (1-3), 6.1 (4, the continuum of ion 5) and the bare nucleus (5)
        lv6(level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1, 2, 0, 1, level, 6], f32[E, g, 1, 54.4])
        two = Radix.level_table([lv(1, 0.0, 2), lv(2, 10.2, 8), lv(3, 12.0, 8), lv(4, 13.6, 1), lv6(1, 0.0, 1), lv6(2, 54.4, 1)];
            masses=Dict(5 => 4.0, 6 => 4.0))
        ions = [Radix.Ion("he_i", Int32(1), Int32(2), Int32(5), 24.6f0), Radix.Ion("he_ii", Int32(2), Int32(2), Int32(6), 54.4f0)]
        @test Radix.element_ions(vcat(reverse(ions), [Radix.Atom("he", Int32(2), Int32(2), 0.1f0, 4.0f0)]), 2) == [5, 6]
        t6(i, f, a, b; type=3) = ToyRate(Int8(type), Int32(6), (; init=i, final=f, frate=a, irate=b))
        rates = [[toy(1, 4, 2.0, 1.0; type=7), toy(1, 2, 3.0, 5.0)], [t6(1, 2, 7.0, 11.0; type=7)]]
        A = Radix.element_matrix(rates, two, [5, 6], cell)
        @test size(A) == (5, 5)
        @test A[4, 1] == 2.0 && A[1, 4] == 1.0 && A[2, 1] == 3.0 && A[1, 2] == 5.0     # ionization of ion 5 goes to ion 6
        @test A[5, 4] == 7.0 && A[4, 5] == 11.0                                         # and that of ion 6 to the nucleus
        # a photoionization that leaves ion 6 in its level 1 has `final` beyond the continuum of ion 5
        A = Radix.element_matrix([[toy(1, 5, 2.0, 1.0; type=7)], Radix.AbstractRate[]], two, [5, 6], cell)
        @test A[5, 1] == 2.0 && A[1, 5] == 1.0
        x = Radix.level_populations(Radix.element_matrix(rates, two, [5, 6], cell))
        @test sum(x) ≈ 1 && all(>=(0), x)
        f = Radix.ion_fractions(x, two, [5, 6])
        @test length(f) == 3 && sum(f) ≈ 1 && f[1] ≈ sum(x[1:3]) && f[2] ≈ x[4] && f[3] ≈ x[5]
        # a level without any rate has no population, and the rest is still solved
        x = Radix.level_populations([-3.0 5.0 0.0; 3.0 -5.0 0.0; 0.0 0.0 0.0][[1, 3, 2], [1, 3, 2]])
        @test x ≈ [5/8, 0, 3/8]
    end

    @testset "the populations" begin
        # two levels: x₂/x₁ is the ratio of the up and the down rates; the populations add up to 1
        x = Radix.level_populations([-3.0 5.0; 3.0 -5.0])
        @test x ≈ [5/8, 3/8]
        # three levels in a chain: detailed balance gives the product of the ratios
        A = Radix.rate_matrix([toy(1, 2, 2.0, 4.0), toy(2, 3, 3.0, 6.0)], levels, 5, cell)
        x = Radix.level_populations(A[1:3, 1:3])
        @test sum(x) ≈ 1 && x[2]/x[1] ≈ 0.5 && x[3]/x[2] ≈ 0.5
        # in LTE
        xl = Radix.lte_populations(levels, 5, cell)
        @test sum(xl) ≈ 1
        @test xl[2]/xl[1] ≈ 4*exp(-Float64(f32(10.2))/Radix.constants().kT_eV) rtol=1e-12
    end
  end
end

# the H I ion and the elements of the database: the matrix reaches LTE at a high density and reproduces the neutral
# fraction of the XSTAR reference run (test/reference/xstar_pow_xi2)
function element_balance_tests(db)
    @testset "Element level balance" begin
        levels = Radix.level_table(db)
        rates = Radix.ion_rates(db, 1)
        @test length(rates) == 550     # (the 33 levels of the ion are not rates)
        @testset "LTE at a high density" begin
            cell = Radix.Cell(1.0, 0.0, 1e22, 1e22)
            x = Radix.level_populations(Radix.rate_matrix(rates, levels, 1, cell))
            xl = Radix.lte_populations(levels, 1, cell)
            # (the constants of XSTAR's rates are rounded differently from those of levwk: 0.5 % here, 10⁻⁵ with CODATA)
            @test x[end] ≈ xl[end] rtol=1e-2
            @test x[2]/x[1] ≈ xl[2]/xl[1] rtol=2e-3
        end
        @testset "the ion fractions of the XSTAR reference run" begin
            # the incident spectrum of the run in 10³⁸ erg/s per erg. Its first zone is at 10¹³ cm (the second row of the
            # table; the first is not converged) and the next, of 10¹³ cm, is evaluated at its midpoint (the table
            # gives its outer edge). The ions with a fraction below `floor` are left out: XSTAR sets those to zero.
            dir = joinpath(@__DIR__, "reference", "xstar_pow_xi2")
            spectrum = fits(joinpath(dir, "xout_cont1.fits"))[3].data
            L38 = 1e38
            E, LE = Float64.(spectrum.energy), Float64.(spectrum.incident)*L38
            abundances = fits(joinpath(dir, "xout_abund1.fits"))[2].data
            roman = ["i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x"]
            floor = 1e-9
            cell = Radix.Cell(100.0, 0.0, 1e4, 1e4)
            for (symbol, Z) in (("h", 1), ("he", 2), ("c", 6), ("n", 7), ("o", 8), ("ne", 10))
                ions = Radix.element_ions(db, Z)
                rates = [Radix.ion_rates(db, ion) for ion in ions]
                for (r, row) in ((1e13, 2), (1.5e13, 3))
                    radiation = Radiation(E, LE/(4π*r^2))
                    x = Radix.level_populations(Radix.element_matrix(rates, levels, ions, cell; radiation))
                    fractions = Radix.ion_fractions(x, levels, ions)
                    for k in eachindex(ions)
                        ref = getproperty(abundances, Symbol(symbol, "_", roman[k]))[row]
                        ref > floor && @test fractions[k] ≈ ref rtol=3e-2
                    end
                end
            end
        end
    end
end

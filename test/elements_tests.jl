# The level balance of one ion: the rate matrix and its populations.

using FITSFiles
using ForwardDiff
using SparseArrays

# a rate between two levels of a toy ion that returns what it is told: (init, final, forward, inverse)
struct ToyRate <: Radix.AbstractRate
    rtype::Int8
    ion::Int32
    r::NamedTuple{(:init, :final, :frate, :irate), Tuple{Int, Int, Float64, Float64}}
end
Radix.rate(c::ToyRate, cell; kw...) = c.r

# a rate between two levels that depends on the temperature: x₂/x₁ = frate/irate = T
struct HotToyRate <: Radix.AbstractRate
    rtype::Int8
    ion::Int32
end
Radix.rate(c::HotToyRate, cell; kw...) = (; init=1, final=2, frate=cell.T^2, irate=cell.T)

function toy_balance_tests()
  @testset "Level balance" begin
    f32 = Float32
    lv(level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1, 2, 0, 1, level, 5], f32[E, g, 1, 13.6])
    levels = Radix.levels([lv(1, 0.0, 2), lv(2, 10.2, 8), lv(3, 12.0, 8), lv(4, 13.6, 1)]; masses=Dict(5 => 1.0))
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
        two = Radix.levels([lv(1, 0.0, 2), lv(2, 10.2, 8), lv(3, 12.0, 8), lv(4, 13.6, 1), lv6(1, 0.0, 1), lv6(2, 54.4, 1)];
            masses=Dict(5 => 4.0, 6 => 4.0))
        ions = [Radix.Ion("he_i", Int32(1), Int32(2), Int32(5), 24.6f0), Radix.Ion("he_ii", Int32(2), Int32(2), Int32(6), 54.4f0)]
        @test Radix.element_ions(vcat(reverse(ions), [Radix.Atom("he", Int32(2), Int32(2), 0.1f0, 4.0f0)]), 2) == [5, 6]
        t6(i, f, a, b; type=3) = ToyRate(Int8(type), Int32(6), (; init=i, final=f, frate=a, irate=b))
        rates = [[toy(1, 4, 2.0, 1.0; type=7), toy(1, 2, 3.0, 5.0)], [t6(1, 2, 7.0, 11.0; type=7)]]
        layout = Radix.Elements(rates, two, [5, 6])
        @test layout.N == 5 && layout.offset == [0, 3] && layout.nlev == [4, 2] && size(layout) == (5, 5)
        A = Radix.element_matrix(layout, cell)
        @test size(A) == (5, 5)
        @test A[4, 1] == 2.0 && A[1, 4] == 1.0 && A[2, 1] == 3.0 && A[1, 2] == 5.0     # ionization of ion 5 goes to ion 6
        @test A[5, 4] == 7.0 && A[4, 5] == 11.0                                         # and that of ion 6 to the nucleus
        # a photoionization that leaves ion 6 in its level 1 has `final` beyond the continuum of ion 5
        A = Radix.element_matrix(Radix.Elements([[toy(1, 5, 2.0, 1.0; type=7)], Radix.AbstractRate[]], two, [5, 6]), cell)
        @test A[5, 1] == 2.0 && A[1, 5] == 1.0
        A0 = Radix.element_matrix(layout, cell)
        x = Radix.level_populations(A0)
        @test sum(x) ≈ 1 && all(>=(0), x)
        f = Radix.ion_fractions(x, layout)
        @test length(f) == 3 && sum(f) ≈ 1 && f[1] ≈ sum(x[1:3]) && f[2] ≈ x[4] && f[3] ≈ x[5]
        # the layout is reused: a second cell refills the same matrix, and the matrix can be given
        A = zeros(5, 5)
        @test Radix.element_matrix!(A, layout, Radix.Cell(2.0, 0.0, 1e4, 1e4)) === A && A == A0
        @test_throws DimensionMismatch Radix.element_matrix!(zeros(4, 4), layout, cell)
        # a level without any rate has no population, and the rest is still solved
        x = Radix.level_populations([-3.0 5.0 0.0; 3.0 -5.0 0.0; 0.0 0.0 0.0][[1, 3, 2], [1, 3, 2]])
        @test x ≈ [5/8, 0, 3/8]
    end

    @testset "escape probabilities and radiation" begin
        line = Radix.AtomicLine2(Int32(4), "", Int32[2, 1, 1, 5], f32[1215.67, 0.4162, 6.265e8], levels)
        layout = Radix.Elements([[line]], levels, [5])
        decay(A) = A[1, 2]
        thin = decay(Radix.element_matrix(layout, cell))
        # `escape` as a pair for each record, or as a function of the record
        @test decay(Radix.element_matrix(layout, cell; escape=[(0.1, 0.2)])) ≈ thin*0.3
        @test decay(Radix.element_matrix(layout, cell; escape=coef -> (0.25, 0.25))) ≈ thin*0.5
        @test decay(Radix.element_matrix(layout, cell; escape=[(0.5, 0.5)])) ≈ thin
        # the radiation: the same flux as a mean intensity or as the flux of a point source
        E = Radix.xstar_energy_grid()
        L = 1e30 ./ E
        r = 1e13
        rad = Radix.point_source(E, L, r)
        @test rad.F ≈ L ./ (4π*r^2) && Radix.mean_intensity(rad) ≈ L ./ (16π^2*r^2)
        @test Radix.Radiation(E; J=Radix.mean_intensity(rad)).F ≈ rad.F
        @test rad.E === E                                     # the grid is shared
        excitation(A) = A[2, 1]
        @test excitation(Radix.element_matrix(layout, cell)) == 0
        @test excitation(Radix.element_matrix(layout, cell; radiation=rad)) > 0
        @test excitation(Radix.element_matrix(layout, cell; radiation=Radix.Radiation(E, 2 .* rad.F))) ≈
            2*excitation(Radix.element_matrix(layout, cell; radiation=rad))
    end

    @testset "a grid of cells" begin
        # every cell has its own radiation and line trapping: the populations of the cell are those of the single-cell solve
        E = Radix.xstar_energy_grid()
        line = Radix.AtomicLine2(Int32(4), "", Int32[2, 1, 1, 5], f32[1215.67, 0.4162, 6.265e8], levels)
        layout = Radix.Elements([[line, toy(1, 3, 1e-2, 5e-2), toy(3, 4, 1e-3, 2e-3)]], levels, [5])
        cells = [Radix.Cell(T, 0.0, 1e4, 1e4) for T in (0.5, 1.0, 2.0), _ in 1:2]
        spectra = [Radix.point_source(E, 1e30 ./ E, 1e13*k) for k in 1:length(cells)]
        trapping = i -> [(0.5/i, 0.5/i), (0.5, 0.5), (0.5, 0.5)]
        x = Radix.element_populations(layout, cells; radiation=i -> spectra[i], escape=trapping)
        @test size(x) == (layout.N, 6)
        for i in eachindex(cells)
            xi = Radix.level_populations(Radix.element_matrix(layout, cells[i]; radiation=spectra[i], escape=trapping(i)))
            @test x[:, i] == xi
        end
        # a fixed value is shared by all the cells
        y = Radix.element_populations(layout, cells; radiation=spectra[1])
        @test y[:, 1] == Radix.level_populations(Radix.element_matrix(layout, cells[1]; radiation=spectra[1]))
    end

    @testset "sparse matrices" begin
        rates = [toy(1, 2, 3.0, 5.0), toy(2, 3, 1.0, 2.0), toy(3, 4, 7.0, 11.0; type=7)]
        layout = Radix.Elements([rates], levels, [5])
        Ad = Radix.element_matrix(layout, cell)
        As = Radix.element_matrix(layout, cell; sparse=true)
        @test As isa SparseMatrixCSC && Matrix(As) == Ad
        @test length(nonzeros(As)) == length(layout.rowval) && all(i -> As[i, i] != 0, 1:4)    # the pattern holds the diagonal
        @test Radix.level_populations(As) ≈ Radix.level_populations(Ad) rtol=1e-12
        # the matrix is refilled in place, and only a matrix with the pattern of the layout is accepted
        @test Radix.element_matrix!(As, layout, Radix.Cell(2.0, 0.0, 1e4, 1e4)) === As
        @test_throws DimensionMismatch Radix.element_matrix!(sparse(zeros(4, 4)), layout, cell)
        # a level without any rate has no population
        x = Radix.level_populations(sparse([-3.0 5.0 0.0; 3.0 -5.0 0.0; 0.0 0.0 0.0][[1, 3, 2], [1, 3, 2]]))
        @test x ≈ [5/8, 0, 3/8]
        # the grid solver gives the same populations with either matrix
        cells = [Radix.Cell(T, 0.0, 1e4, 1e4) for T in (0.5, 1.0, 2.0)]
        @test Radix.element_populations(layout, cells; sparse=true) ≈ Radix.element_populations(layout, cells; sparse=false) rtol=1e-12
        # automatic: sparse from `sparse_from` unknowns on
        @test Radix.sparse_from > layout.N
    end

    @testset "automatic differentiation of the populations" begin
        # x₂ = T/(1 + T): dx₂/dT = 1/(1 + T)²
        pair = Radix.levels([lv(1, 0.0, 2), lv(2, 13.6, 1)]; masses=Dict(5 => 1.0))
        layout = Radix.Elements([[HotToyRate(Int8(3), Int32(5))]], pair, [5])
        function x2(T, sparse)
            A = Radix.element_matrix(layout, Radix.Cell(T, zero(T), one(T), one(T)); sparse)
            Radix.level_populations(A)[2]
        end
        for sparse in (false, true), T in (0.5, 2.0)
            @test ForwardDiff.derivative(t -> x2(t, sparse), T) ≈ 1/(1 + T)^2 rtol=1e-10
        end
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

# the elements of the database: the matrix reaches LTE at a high density and reproduces the neutral
# fraction of the XSTAR reference run (test/reference/xstar_pow_xi2)
function element_balance_tests(db)
    @testset "Element level balance" begin
        levels = Radix.levels(db)
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
        @testset "sparse and dense matrices agree" begin
            cell = Radix.Cell(100.0, 0.0, 1e4, 1e4)
            for Z in (2, 8)
                layout = Radix.Elements(db, levels, Z)
                dense = Radix.element_matrix(layout, cell)
                sparse = Radix.element_matrix(layout, cell; sparse=true)
                @test Matrix(sparse) == dense
                Z == 8 && @test length(nonzeros(sparse)) < length(dense) ÷ 10
                # (the iterative refinement gives the tiny ion fractions too)
                @test Radix.ion_fractions(Radix.level_populations(sparse), layout) ≈
                      Radix.ion_fractions(Radix.level_populations(dense), layout) rtol=1e-6
            end
        end
        @testset "the layout agrees with the rates" begin
            # the matrix entries of a record are those of its levels at any temperature (index=true), including the records
            # that exist only above a temperature
            cell = Radix.Cell(100.0, 0.0, 1e4, 1e4)
            layout = Radix.Elements(db, levels, 2)
            missing_pairs = 0
            for (j, coef) in enumerate(layout.rates)
                r = Radix.rate(coef, cell)
                (r.init > 0 && r.final > 0 && r.init != r.final) || continue
                k, offset = layout.ion[j], layout.offset[layout.ion[j]]
                (coef.rtype in Radix.unordered_types || max(r.init, r.final) <= layout.nlev[k]) || continue
                layout.lo[j] == 0 && (missing_pairs += 1)
            end
            @test missing_pairs == 0
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
            roman = ["i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x", "xi", "xii", "xiii", "xiv", "xv", "xvi", "xvii", "xviii",
                     "xix", "xx", "xxi", "xxii", "xxiii", "xxiv", "xxv", "xxvi"]
            floor, loose_below = 1e-9, 1e-5
            cell = Radix.Cell(100.0, 0.0, 1e4, 1e4)
            for (symbol, Z) in (("h", 1), ("he", 2), ("c", 6), ("n", 7), ("o", 8), ("ne", 10), ("fe", 26))
                layout = Radix.Elements(db, levels, Z)
                for (r, row) in ((1e13, 2), (1.5e13, 3))
                    radiation = Radix.point_source(E, LE, r)
                    x = Radix.level_populations(Radix.element_matrix(layout, cell; radiation, sparse=layout.N >= Radix.sparse_from))
                    fractions = Radix.ion_fractions(x, layout)
                    for k in eachindex(layout.ions)
                        ref = getproperty(abundances, Symbol(symbol, "_", roman[k]))[row]
                        # (XSTAR converges the fractions to its `critf` = 10⁻⁷ of the element only, so the small ones are looser)
                ref > floor && @test fractions[k] ≈ ref rtol=(ref > loose_below ? 3e-2 : 1e-1)
                    end
                end
            end
        end
    end
end

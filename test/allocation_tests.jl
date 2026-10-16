# The rates allocate nothing once they are compiled: their tables are read as they are (converted on access), their small temporaries are tuples and
# the work arrays of the integrals over the grid belong to the task. The records and cells of the ucalc fixtures (test/reference/ucalc).

allocated(coef, cell, kw) = @allocated Radix.rate(coef, cell; kw...)

@testset "The rates do not allocate" begin
    @testset "views and work arrays" begin
        v = Float32[1.5, 2.5, 3.5]
        w = Radix.as_float64(v)
        @test w == Float64.(v) && eltype(w) == Float64 && size(w) == (3,)
        @test Radix.as_float64(Float64[1, 2]) isa Vector{Float64}                       # (it is itself)
        m = Radix.as_type(Float32, rand(2, 3))
        @test size(m) == (2, 3) && eltype(m) == Float32 && m[2, 3] isa Float32
        a = Radix.scratch(:test_a, 5)
        a .= 1:5
        @test length(a) == 5 && Radix.scratch(:test_a, 5) === a                         # (the same array)
        @test length(Radix.scratch(:test_a, 3)) == 3 && length(Radix.scratch(:test_b, 2, Bool)) == 2
        @test Radix.scaled_table([1.0, 2.0, 3.0], 2.0, :test_c; floor=2.5, cap=5.0) == [2.5, 4.0, 5.0]
        @test Radix.scaled_table([1.0, 2.0], 2.0 + 0im, :test_c) isa Vector{ComplexF64}   # (other numbers: a new array)
        c = Radix.Capped(Float64[1, 9, 3], 2, 8.0)
        @test c == [1.0, 8.0, 3.0] && Radix.Capped(Float64[1, 5, 3], 2, 8.0)[2] == 5.0
    end

    @testset "the rates of the fixtures" begin
        for name in ("type01", "type02", "type07", "type09", "type10", "type22", "type30", "type38", "type39", "type51", "type54", "type56", "type57", "type60", "type62",
                "type63", "type65", "type66", "type67", "type68", "type69", "type71", "type72", "type73", "type75", "type76", "type77", "type79", "type81", "type86", "type89", "type92", "type95", "type96", "type98", "type102")
            kw0 = name in ("type02", "type09", "type10") ? (; nlev=5) : NamedTuple()
            worst = 0
            for (c, _) in Iterators.take(ucalc_cases(name), 20)
                coef, cell, lv = ucalc_inputs(c)
                kw = name == "type10" ? kw0 : kw0
                cell2 = name == "type10" ? Radix.Cell(cell.T, cell.ntot - c.cond[5], cell.nₑ, cell.ntot) : cell
                allocated(coef, cell2, kw)
                worst = max(worst, allocated(coef, cell2, kw))
            end
            @test worst == 0
        end
    end

    @testset "the rates on the grid, with the radiation and the opacity arrays of the caller" begin
        for name in ("type49", "type53", "type59", "type64", "type70", "type85", "type88", "type99", "type74", "type50", "type91")
            worst = 0
            for (c, _) in Iterators.take(ucalc_cases(name), 12)
                coef, cell, lv = ucalc_inputs(c)
                rad = ucalc_radiation(c)
                op = Radix.Opacity(9999)
                kw = if name in ("type50", "type91")
                    (; radiation=rad, vturb=c.cond[6], pesc=c.cond[11] + c.cond[12], cfrac=c.cond[8])
                elseif name == "type74"
                    (; radiation=rad)
                else
                    (; radiation=rad, opacity=op, ptmp=(c.cond[11], c.cond[12]), abund=(c.cond[9], c.cond[10]), lfast=Int(c.rad[8]), (name == "type70" ? (; neutral=true) : (;))...)
                end
                allocated(coef, cell, kw)
                worst = max(worst, allocated(coef, cell, kw))
            end
            @test worst == 0
        end
    end

    @testset "the exact binning of the edges and the lines" begin
        # (the constants of ucalc, the comparisons with XSTAR above, do not have `exact_binning`)
        Radix.constants(Radix.Constants()) do
            for name in ("type49", "type53", "type59")
                worst = 0
                for (c, _) in Iterators.take(ucalc_cases(name), 12)
                    coef, cell, lv = ucalc_inputs(c)
                    rad = ucalc_radiation(c)
                    op = Radix.Opacity(9999)
                    kw = (; radiation=rad, opacity=op, ptmp=(c.cond[11], c.cond[12]), abund=(c.cond[9], c.cond[10]), lfast=Int(c.rad[8]))
                    allocated(coef, cell, kw)
                    worst = max(worst, allocated(coef, cell, kw))
                    kw = (; kw..., integrals=false)
                    name == "type59" || (allocated(coef, cell, kw); worst = max(worst, allocated(coef, cell, kw)))
                end
                @test worst == 0
            end
            E = Radix.xstar_energy_grid(10000)
            rad = Radix.Radiation(E, zeros(10000))
            op = Radix.Opacity(10000)
            line!(op, rad, a) = Radix.add_line!(op, rad, 1.0, 0.0, 0.0, 12.0, 1.0, 10.0, 16.0, a)
            for a in (1e-8, 1e-3, 0.5)
                line!(op, rad, a)
                @test @allocated(line!(op, rad, a)) == 0
            end
        end
    end
end

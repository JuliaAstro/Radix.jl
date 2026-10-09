# The opacity of the lines in the continuum bins: the Voigt function and linopac against the real XSTAR (test/reference/ucalc/drvlinopac.f90),
# the ranking of the lines in a bin and the lines in the continuum opacity.

function lines_tests()
    @testset "Lines in the continuum bins" begin
        # the real linopac and voigte of XSTAR: the input is the grid, the lines (optpp, rcem1, rcem2, elin, vturbi, t, aatmp, delea, lfast)
        # and the pairs (v, a) of the Voigt function; the output the bins that changed, and the values of voigte
        data = joinpath(@__DIR__, "reference", "ucalc", "data")
        tokens = split(read(joinpath(data, "linopac.in"), String))
        n = parse(Int, tokens[1])
        E = parse.(Float64, tokens[2:n + 1])
        nlines = parse(Int, tokens[n + 2])
        columns = 9
        lines = [parse.(Float64, tokens[n + 3 + columns*(i - 1):n + 2 + columns*i]) for i in 1:nlines]
        reference = [Tuple{Int, Float64, Float64, Float64}[] for _ in 1:nlines]
        voigts = NTuple{3, Float64}[]
        current = 0
        for text in eachline(joinpath(data, "linopac.out"))
            fields = split(text)
            if fields[1] == "line"
                current = parse(Int, fields[2])
            elseif fields[1] == "voigt"
                push!(voigts, Tuple(parse.(Float64, fields[2:4])))
            else
                push!(reference[current], (parse(Int, fields[1]), parse.(Float64, fields[2:4])...))
            end
        end

        @testset "voigt" begin
            for (v, a, h) in voigts
                @test Radix.voigt(v, a) == h                                  # the same single-precision constants: identical
            end
            @test Radix.voigt(-0.7, 0.05) == Radix.voigt(0.7, 0.05)
            @test Radix.voigt(0.0, 0.0) == 1 && Radix.voigt(200.0, 0.0) == 0
            @test Radix.voigt(0.0, 1e-3) < 1 && Radix.voigt(0.0, 1e-3) > 0.99
            # the profile is normalized: the integral over v is √π (the damping wings are slow: a long range)
            for a in (0.0, 0.01, 0.1)
                @test sum(Radix.voigt(v, a) for v in -400:0.01:400)*0.01 ≈ sqrt(pi) rtol=2e-3
            end
        end

        @testset "linopac" begin
            Radix.with_constants(Radix.ucalc_constants()) do
                rad = Radix.Radiation(E, zeros(n))
                for (i, (optpp, rcem1, rcem2, elin, vturb, T, mass, delea, lfast)) in enumerate(lines)
                    opacity = Radix.Opacity(n)
                    Radix.add_line!(opacity, rad, optpp, rcem1, rcem2, elin, vturb, T, mass, delea; lfast=Int(lfast))
                    changed = [r[1] for r in reference[i]]
                    @test findall(!=(0), vec(sum(abs, opacity.emissivity; dims=1)) .+ opacity.total) == changed      # the same bins
                    @test all(r -> isapprox(opacity.total[r[1]], r[2]; rtol=1e-6), reference[i])
                    @test all(r -> isapprox(opacity.emissivity[1, r[1]], r[3]; rtol=1e-6) && isapprox(opacity.emissivity[2, r[1]], r[4]; rtol=1e-6), reference[i])
                end
            end
        end

        @testset "the total opacity of a line" begin
            # the bins add up to the integral of the profile, `optpp × width`, when the line is well inside the grid
            Radix.with_constants(Radix.ucalc_constants()) do
                rad = Radix.Radiation(E, zeros(n))
                optpp, elin, T, mass = 1e-3, 12398.4016/1500.0, 1.0, 56.0
                for delea in (1e-8, 1e-3)
                    opacity = Radix.Opacity(n)
                    Radix.add_line!(opacity, rad, optpp, 0.0, 0.0, elin, 100.0, T, mass, delea)
                    K = Radix.constants()
                    e0 = K.hc_eVÅ_single/elin
                    width = sqrt((e0*sqrt(T/mass)*K.thermal_speed/Radix.cm_per_km/(K.light_speed/Radix.cm_per_km))^2 + (e0*100.0/(K.light_speed/Radix.cm_per_km))^2)
                    area = sum(opacity.total[k]*(E[k + 1] - E[k]) for k in 1:n - 1)
                    @test area ≈ optpp*width rtol=0.05
                end
                # a single bin (lfast > 2): the emissivities too
                opacity = Radix.Opacity(n)
                Radix.add_line!(opacity, rad, 1.0, 2.0, 3.0, 12398.4016/1500.0, 100.0, 1.0, 56.0, 1e-8; lfast=3)
                @test count(!=(0), opacity.total) == 1 && count(!=(0), opacity.emissivity[1, :]) == 1
                @test opacity.emissivity[2, argmax(opacity.total)] ≈ 1.5*opacity.emissivity[1, argmax(opacity.total)]
                # outside the range of the wavelengths, and below the grid, nothing is added
                none = Radix.Opacity(n)
                Radix.add_line!(none, rad, 1.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 0.0)
                Radix.add_line!(none, rad, 1.0, 0.0, 0.0, 2e8, 1.0, 1.0, 1.0, 0.0)
                Radix.add_line!(none, rad, 1.0, 0.0, 0.0, 12398.4016/E[1]*2, 1.0, 1.0, 1.0, 0.0)
                @test all(==(0), none.total)
            end
        end

        @testset "the natural width of a line" begin
            f32 = Float32
            lv(level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1, 2, 0, 1, level, 5], f32[E, g, 1, 13.6])
            levels = Radix.levels([lv(1, 0.0, 2), lv(2, 10.2, 8), lv(3, 13.6, 1)]; masses=Dict(5 => 1.0))
            line = Radix.AtomicLine2(Int32(4), "", Int32[2, 1, 1, 5], f32[1215.67, 0.4162, 6.265e8], levels)
            other = Radix.AtomicLine2(Int32(4), "", Int32[3, 1, 1, 5], f32[912.0, 0.1, 1e8], levels)
            vacancy = Radix.IronKAuger(Int32(41), "", Int32[1, 2, 26, 4, 5], f32[10.2, 1e14, 5e13, 1e9], levels)     # the level 2 is a K-vacancy level
            second = Radix.IronKAuger(Int32(41), "", Int32[1, 2, 26, 4, 5], f32[10.2, 2e14, 9e13, 1e9], levels)      # (the first record of a level is the one that deleafnd takes)
            mixture = Radix.Mixture(levels, [26], [1.0], [Radix.Elements([[line, other, vacancy, second]], levels, [5])])
            widths = Radix.auger_widths(mixture)
            @test widths == Dict((5, 2) => Float64(f32(5e13)))
            h = Float64(4.136f-15)
            @test Radix.line_width(line) == Float64(f32(6.265e8))*h                              # the Einstein A without the table
            @test Radix.line_width(line, Dict{Tuple{Int, Int}, Float64}()) == Float64(f32(6.265e8))*h
            @test Radix.line_width(line, widths) == Float64(f32(5e13))*h                         # the upper level of the line is the vacancy
            @test Radix.line_width(other, widths) == Float64(f32(1e8))*h                         # another line keeps A
            @test Radix.line_width(vacancy) === nothing                                          # (only the lines are put in the bins)
            # and the profile of the line follows: the damping parameter is 10⁷ times larger, the peak of the bins lower and the wings higher
            E = Radix.xstar_energy_grid(999)
            rad = Radix.Radiation(E, zeros(length(E)))
            narrow, broad = Radix.Opacity(length(E)), Radix.Opacity(length(E))
            Radix.add_line!(narrow, rad, line, 1e-3, (0.0, 0.0), 1.0)
            Radix.add_line!(broad, rad, line, 1e-3, (0.0, 0.0), 1.0; widths)
            @test maximum(broad.total) < maximum(narrow.total)
            k = argmax(narrow.total)
            @test broad.total[k + 5] > narrow.total[k + 5]
        end

        @testset "the ranking of a bin" begin
            list = Tuple{Float64, Int, Int}[]
            for (strength, j) in ((5.0, 1), (7.0, 2), (6.0, 3), (1.0, 4))
                Radix.rank!(list, (strength, 1, j), 4)
            end
            @test [e[3] for e in list] == [2, 3, 1]                  # the last place (nrank) is never filled by an entry that ends there
            Radix.rank!(list, (9.0, 1, 5), 4)
            @test [e[3] for e in list] == [5, 2, 3, 1]               # but it is by one that pushes the others down
            Radix.rank!(list, (4.0, 1, 6), 4)
            @test [e[3] for e in list] == [5, 2, 3, 1]               # and an entry whose place is the last is not stored
            Radix.rank!(list, (6.5, 1, 7), 4)
            @test [e[3] for e in list] == [5, 2, 7, 3]
            Radix.rank!(list, (6.5, 1, 8), 4)
            @test [e[3] for e in list] == [5, 2, 7, 3]               # (6.5 again: it ranks below the other 6.5 and above the last, whose place is the last: not stored)
        end
    end
end

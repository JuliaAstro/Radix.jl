# The level data of a rate is stored in its coefficients (a LevelTable in the field `levels`), so that
# `rate(coef, cell)` needs nothing else.

@testset "Level data in the coefficients" begin
    f32 = Float32
    lv(ion, level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1, 2, 0, 1, level, ion], f32[E, g, 1, 13.6])
    atom = Radix.Atom(Int32(13), "iron", Int32[26, 26], f32[2.5e-5, 55.85])
    ionrec = Radix.Ion(Int32(14), "fe_i", Int32[1, 26, 7], f32[7.9])
    records = [lv(7, 1, 0.0, 2), lv(7, 2, 10.2, 8), lv(7, 3, 13.6, 1), lv(8, 1, 0.0, 1), lv(8, 2, 5.0, 3), lv(8, 3, 9.0, 1), atom, ionrec]
    table = Radix.level_table(records)

    @testset "the table" begin
        @test table isa LevelTable
        @test length(table) == 6 && haskey(table, (7, 2)) && table[(7, 2)].g == 8     # a dictionary of levels
        @test get(table, (7, 9), nothing) === nothing
        @test nlevels(table, 7) == 3 && nlevels(table, 8) == 3 && nlevels(table, 99) == 0
        @test level_counts(table) == Dict(7 => 3, 8 => 3)
        @test atomic_mass(table, 7) == Float64(f32(55.85))                          # from the Atom and Ion records
        @test_throws ArgumentError atomic_mass(table, 8)                            # ion 8 has no element record
        @test atomic_mass(Radix.level_table(records; masses=Dict(8 => 16.0)), 8) == 16.0
    end

    @testset "attaching" begin
        c = Radix.CollisionFe19(Int8(3), "x", Radix.Transition(Int32(1), Int32(2)), Int32(26), Int32(7), f32(0.5))
        @test c.levels === nothing
        a = Radix.attach_levels(c, table)
        @test a.levels === table && typeof(a).name.wrapper === Radix.CollisionFe19
        @test a.Υ == c.Υ && a.transition == c.transition                             # everything else is kept
        # a rate without level data says what to do
        cell = Radix.Cell(1.0, 1e3, 2e3, 1e4)
        err = try Radix.rate(c, cell) catch e e end
        @test err isa ArgumentError && occursin("attach_levels", err.msg)
        @test Radix.rate(a, cell).frate > 0 && a(cell) == Radix.rate(a, cell)       # the functor needs only the cell
        # rates that need no level data are unchanged by attaching
        rr = Radix.RadRecomb(Int32(1), "rr", Int32[1], f32[1e-12, 0.7])
        @test Radix.attach_levels(rr, table) === rr
        # the number of levels and the mass come from the table
        line = Radix.attach_levels(Radix.AtomicLine2(Int32(4), "x", Int32[2, 1, 1, 7], f32[1215.67, 0.4162, 6.265e8]), table)
        @test Radix.rate(line, cell).init == 2
        @test Radix.rate(line, cell; mass=1.0).frate == Radix.rate(line, cell).frate        # (the mass only enters the opacity)
        @test Radix.rate(line, cell).opacity != Radix.rate(line, cell; mass=1.0).opacity
        nomass = Radix.attach_levels(Radix.AtomicLine2(Int32(4), "x", Int32[2, 1, 1, 8], f32[1215.67, 0.4162, 6.265e8]), table)
        @test_throws ArgumentError Radix.rate(nomass, cell)
        @test Radix.rate(nomass, cell; mass=16.0).frate > 0
    end

    @testset "photoionization without a spectrum" begin
        # by default the photoionization rates see no radiation: no photoionization, and the gas's own recombination
        cell = Radix.Cell(1.0, 1e3, 2e3, 1e4)
        E = Radix.xstar_energy_grid()
        @test length(E) == 9999 && E[1] == 0.1 && 3.99e5 < E[9800] < 4.01e5 && 1.0e6 < E[end] < 1.01e6
        @test Radix.NO_RADIATION.E == E && all(iszero, Radix.NO_RADIATION.F)
        p3 = Radix.attach_levels(Radix.ParPhotoIonize3(Int32(59), "p", Int32[1, 1, 0, 1, 7, 1, 7], f32[13.6, 10.0, 1.0e3, 1.0, 3.0, 1.0]), table)
        @test Radix.rate(p3, cell).frate == 0
    end
end

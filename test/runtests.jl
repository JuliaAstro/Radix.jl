using Radix
using Test

# Rate formulas below are transcribed independently from XSTAR's ucalc()
# (ftools/xstar/utils/xstarsub.f) so Radix is checked against the Fortran,
# not against itself. T is in units of 1e4 K, as in XSTAR.

@testset "Radix.jl" begin
    f32 = Float32
    T, nh, ne = 1.5, 1e3, 2e3
    cell = Radix.Cell(T, nh, ne)

    @testset "rate interface" begin
        coefs = [
            Radix.RadRecomb(Int32(1), "rr", Int32[1], f32[1e-12, 0.7]),
            Radix.Autoionize(Int32(3), "ai", Int32[1], f32[1e-3, 1.0]),
            Radix.AtomicLevel(Int32(6), "lv", Int32[1,1,0,1,1,1], f32[0,2,1,13.6]),
            Radix.TotRadRecombH(Int32(30), "rrH", Int32[2, 3], f32[]),
        ]
        for c in coefs
            r = Radix.rate(c, cell)
            @test keys(r) == (:init, :final, :frate, :irate)
            @test Radix.rate(c, cell; index=true).frate == 0
        end
    end

    @testset "Transition and Parent" begin
        c = Radix.CollisionHlike1(Int32(60), "c", Int32[3,5,1,2], f32[1,2,3])
        @test c.transition == Radix.Transition(Int32(3), Int32(5))
        @test c.ion == 2
        p = Radix.ParPhotoIonize1(Int32(7), "p", Int32[2,1,3,26,81,22,1,21], f32[1,2])
        @test p.parent == Radix.Parent(Int32(22), Int32(81))   # ion, level
        @test (p.level, p.ion) == (1, 21)
    end

    @testset "functor convenience" begin
        c = Radix.RadRecomb(Int32(1), "rr", Int32[1], f32[1e-12, 0.7])
        @test c(cell) == Radix.rate(c, cell)
        @test c(cell; index=true) == Radix.rate(c, cell; index=true)
        cells = [Radix.Cell(T, nh, ne) for T in (0.5, 1.0, 2.0)]
        @test [r.frate for r in c.(cells)] == [Radix.rate(c, x).frate for x in cells]
    end

    @testset "unported rates error clearly" begin
        c = Radix.CollisionHlike1(Int32(60), "c", Int32[1,2,1,0,0,0,0,0], f32[1,2,3])
        @test_throws ErrorException Radix.rate(c, cell)
    end

    @testset "XSTAR reference formulas" begin
        # ucalc label 1: rrrt = arad/t**eta; ans1 = rrrt*xnx
        a, η = 1e-12, 0.7
        c = Radix.RadRecomb(Int32(1), "rr", Int32[1], f32[a, η])
        r = Radix.rate(c, cell)
        @test r.frate ≈ ne*f32(a)/T^f32(η)
        @test (r.init, r.final) == (1, 0)

        # ucalc label 3: airt = cai*expo(-eai/(t*0.861707))/sqrt(t); ans1 = airt*xnx
        cai, eai = 1e-3, 1.0
        c = Radix.Autoionize(Int32(3), "ai", Int32[1], f32[cai, eai])
        r = Radix.rate(c, cell)
        @test r.frate ≈ ne*f32(cai)*exp(-f32(eai)/(T*0.861707))/sqrt(T)
        @test (r.init, r.final) == (1, 1)

        # ucalc label 30 (hydrogenic RR): t6=t/100; beta=Z²/(6.34 t6); vth=3.10782e7 sqrt(t);
        # phi from the two fits joined with a fudge factor; ans1 = 2*2.105e-22*vth*beta*phi*xnx
        Zc = 3
        c = Radix.TotRadRecombH(Int32(30), "rrH", Int32[Zc, 3], f32[])
        t6 = T/100; β = Zc^2/(6.34*t6); vth = 3.10782e7*sqrt(T)
        ypow = min(1.0, 0.06376/β^2)
        fudge = 0.9*(1-ypow) + (1/1.5)*ypow
        ϕ1 = (1.735 + log(β) + 1/6/β)*fudge/2
        ϕ2 = β*(-1.202*log(β) - 0.298)
        ϕ = β < 0.2525 ? ϕ2 : ϕ1
        @test Radix.rate(c, cell).frate ≈ ne*2*2.105e-22*vth*β*ϕ

        # ucalc label 2: rate = a*expo(log(t)*b)*(1+c*expo(d*t))*1e-9; ans1 = rate*xh0
        A, B, C, D = 2.0, 0.5, 0.3, -0.2
        c = Radix.ChargeExH0(Int32(2), "cx", Int32[1], f32[A, B, C, D, 0, 0, 0])
        r = Radix.rate(c, cell)
        @test r.frate ≈ nh*f32(A)*exp(log(T)*f32(B))*(1 + f32(C)*exp(f32(D)*T))*1e-9
    end

    # The XSTAR output files kept as the target for the solver; check they are
    # present and readable (see test/reference/xstar_pow_xi2/README.md).
    @testset "XSTAR reference outputs" begin
        using FITSFiles
        dir = joinpath(@__DIR__, "reference", "xstar_pow_xi2")
        # file => (number of HDUs, [HDU index => columns expected there])
        expected = Dict(
            "xout_abund1.fits" => (5, [2 => [:radius, :temperature, :h_i, :fe_xxvi],
                                       5 => [:radius, :hydrogen, :total]]),
            "xout_cont1.fits"  => (3, [3 => [:energy, :incident, :transmitted]]),
            "xout_lines1.fits" => (3, [3 => [:ion, :wavelength, :emit_outward]]),
            "xout_rrc1.fits"   => (3, [3 => [:ion, :level, :energy]]),
            "xout_spect1.fits" => (3, [3 => [:energy, :incident, :transmitted]]),
        )
        @test isfile(joinpath(dir, "README.md"))
        @test isfile(joinpath(dir, "xstar.par"))
        @test isfile(joinpath(dir, "xout_step.log"))
        for (name, (nhdu, tables)) in expected
            path = joinpath(dir, name)
            @test isfile(path)
            hdus = fits(path)
            @test length(hdus) == nhdu
            for (i, cols) in tables
                data = hdus[i].data
                for c in cols
                    @test haskey(data, c)
                    @test length(data[c]) > 0
                end
                # all columns of a table have the same number of rows
                @test length(unique(length(v) for v in values(data))) == 1
            end
        end
    end

    # Parse the full XSTAR atomic database and compare record counts per type.
    # Slow (~35 s, 870 MB): opt in with RADIX_ATDB=/path/to/atdb.fits
    atdb = get(ENV, "RADIX_ATDB", "")
    if isfile(atdb)
        @testset "atdb.fits parsing" begin
            expected = Dict{String,Int}()
            for line in eachline(joinpath(@__DIR__, "reference", "atdb_counts.txt"))
                startswith(line, "#") && continue
                k, v = split(line)
                expected[k] = parse(Int, v)
            end
            db = open(Radix.load, atdb)
            got = Dict{String,Int}()
            for r in db
                k = String(nameof(typeof(r)))
                got[k] = get(got, k, 0) + 1
            end
            @test got == expected
        end
    else
        @info "Skipping atdb.fits parsing test (set RADIX_ATDB to enable)"
    end
end

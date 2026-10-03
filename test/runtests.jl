using Radix
using Test

include("ucalc_tests.jl")
include("forwarddiff_tests.jl")

# Rate formulas below are transcribed independently from XSTAR's ucalc()
# (ftools/xstar/utils/xstarsub.f) so Radix is checked against the Fortran,
# not against itself. T is in units of 1e4 K, as in XSTAR.

@testset "Radix.jl" begin
    f32 = Float32
    T, nh, ne = 1.5, 1e3, 2e3
    cell = Radix.Cell(T, nh, ne, 1e4)

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

    @testset "AtomicLine2" begin
        # numbers are checked against ucalc in ucalc_tests.jl; here the edge cases
        lv(level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1,2,0,1,level,7],
            f32[E, g, 1, 13.6])
        levels = Radix.level_table([lv(1, 0.0, 2), lv(2, 10.2, 8)])
        mk(i, k; λ=1215.67, A=6.265e8, rt=4) = Radix.AtomicLine2(Int32(rt), "x", Int32[i,k,1,7], f32[λ, 0.4162, A])
        c = mk(2, 1)
        r = Radix.rate(c, cell; levels=levels, mass=1.01)
        @test (r.init, r.final) == (2, 1)                 # upper first, whatever the stored order
        @test Radix.rate(mk(1, 2), cell; levels=levels, mass=1.01) == r
        @test r.irate == 0 && r.ienergy == 0              # no photoexcitation without a spectrum
        # photoexcitation from the spectrum at the line energy, reduced by the covering fraction
        E = [0.1*exp(0.0015*(i - 1)) for i in 1:9999]
        rad = Radix.Radiation(E, 1e16 .* E .^ -1.0)
        rp = Radix.rate(c, cell; levels=levels, mass=1.01, radiation=rad)
        @test rp.irate > 0 && rp.frate == r.frate
        @test rp.ienergy ≈ rp.irate*10.2*1.602176634e-12  rtol=1e-6    # ΔE of the two levels
        @test Radix.rate(c, cell; levels=levels, mass=1.01, radiation=rad, cfrac=0.25).irate ≈ 0.75*rp.irate
        @test Radix.rate(mk(2, 1; λ=1e9), cell; levels=levels, mass=1.01, radiation=rad).irate == 0
        @test r.frate == Float64(f32(6.265e8))            # A times the escape probability (1)
        @test Radix.rate(c, cell; levels=levels, mass=1.01, pesc=0.5).frate ≈ r.frate/2
        # a tiny A is floored at 1e-20 times the density
        @test Radix.rate(mk(2, 1; A=1e-30), cell; levels=levels, mass=1.01).frate ≈ 1e-20*cell.ntot
        # lines without a wavelength give nothing; very long wavelengths have no opacity
        @test Radix.rate(mk(2, 1; λ=0), cell; levels=levels, mass=1.01).frate == 0
        far = Radix.rate(mk(2, 1; λ=1e9), cell; levels=levels, mass=1.01)
        @test far.frate > 0 && far.opacity == 0
        # unknown levels, and the continuum level (index nlev), give nothing
        @test Radix.rate(mk(2, 9), cell; levels=levels, mass=1.01).frate == 0
        @test Radix.rate(c, cell; levels=levels, mass=1.01, nlev=2).frate == 0
        @test Radix.rate(c, cell; index=true, levels=levels, mass=1.01).frate == 0
    end

    @testset "ElectronImpact1" begin
        lv(level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1,2,0,1,level,7],
            f32[E, g, 1, 13.6])
        levels = Radix.level_table([lv(1, 0.0, 2), lv(2, 10.2, 8)])
        Tg, Ug = [3.0, 4.0, 5.0], [0.1, 0.4, 0.2]          # log10 T (K), Υ
        mk(i, k) = Radix.ElectronImpact1(Int8(5), "x", Radix.Transition(Int32(i), Int32(k)),
            Int32(1), Int32(7), f32.(Tg), f32.(Ug))
        # ucalc label 56 (T in 1e4 K, tfnd = log10(1e4 t)):
        #   cij = 8.626e-8 Υ exp(-ΔE/(0.861707 t))/(sqrt(t) gglo), cji = 8.626e-8 Υ/(sqrt(t) ggup)
        Tc = 0.5                                              # log10(5000 K) = 3.699
        Υi = 0.1 + (0.4 - 0.1)*(log10(Tc*1e4) - 3.0)/(4.0 - 3.0)
        c = mk(1, 2)
        r = Radix.rate(c, Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels)
        Δ = Float64(f32(10.2))
        cij = 8.626e-8*Υi*exp(-Δ/(0.861707*Tc))/sqrt(Tc)/2
        cji = 8.626e-8*Υi/sqrt(Tc)/8
        @test (r.init, r.final) == (1, 2)
        @test r.frate ≈ cij*2e3  rtol=1e-6
        @test r.irate ≈ cji*2e3  rtol=1e-6
        # detailed balance: cij/cji = (gup/glo) exp(-ΔE/kT)
        @test r.frate/r.irate ≈ (8/2)*exp(-Δ/(0.861707*Tc))  rtol=1e-6
        # stored upper-first gives the same answer
        @test Radix.rate(mk(2, 1), Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels) == r
        # above the table: last segment extrapolated and clamped at 0
        hi = Radix.rate(c, Radix.Cell(10.0, 1e3, 2e3, 1e4); levels=levels)   # log T = 5.0
        @test hi.irate ≈ 8.626e-8*0.2/sqrt(10.0)/8*2e3  rtol=1e-6
        far = Radix.rate(c, Radix.Cell(1e3, 1e3, 2e3, 1e4); levels=levels)    # log T = 7.0
        @test far.frate == 0 && far.irate == 0               # 0.4 + (0.2-0.4)*3 < 0
        # unknown level: nothing
        @test Radix.rate(mk(1, 9), Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels).frate == 0
    end

    @testset "ElectronCollision" begin
        # the 5-point spline passes through its knots and is exact for linear data
        kn = [0.3, 1.1, 0.7, 2.0, 1.5]
        @test [Radix.splinem(kn, x) for x in (0.0, 0.25, 0.5, 0.75, 1.0)] ≈ kn
        @test all(x -> Radix.splinem(0:0.25:1, x) ≈ x, 0:0.05:1)
        # upsil(): x = log((e+c)/c)/log(e+c) (kinds 1,4) or e/(e+c) (2,3), then the kind's factor
        e(T, ΔE) = abs(T/(1.57888e5*ΔE))
        T0, ΔE0, C0 = 2e5, 0.5, 1.7
        ee = e(T0, ΔE0)
        @test Radix.chianti_upsilon(2, ΔE0, C0, kn, T0) ≈ Radix.splinem(kn, ee/(ee + C0))
        @test Radix.chianti_upsilon(1, ΔE0, C0, kn, T0) ≈
            Radix.splinem(kn, log((ee + C0)/C0)/log(ee + C0))*log(ee + 2.71828)
        @test Radix.chianti_upsilon(3, ΔE0, C0, kn, T0) ≈ Radix.splinem(kn, ee/(ee + C0))/(ee + 1)
        @test Radix.chianti_upsilon(4, ΔE0, C0, kn, T0) ≈
            Radix.splinem(kn, log((ee + C0)/C0)/log(ee + C0))*log(ee + C0)
        @test_throws ArgumentError Radix.chianti_upsilon(5, ΔE0, C0, kn, T0)

        lv(level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1,2,0,1,level,7],
            f32[E, g, 1, 13.6])
        levels = Radix.level_table([lv(1, 0.0, 2), lv(2, 6.8, 8)])
        mk(i, k; ΔE=0.5) = Radix.ElectronCollision(Int8(3), "x", Int32(2),
            Radix.Transition(Int32(i), Int32(k)), Int32(1), Int32(7),
            f32(ΔE), f32(C0), f32.(kn))
        # the numbers are checked against ucalc in ucalc_tests.jl
        Tc = 20.0
        c = mk(1, 2)
        r = Radix.rate(c, Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels)
        @test (r.init, r.final) == (1, 2)
        @test r.frate/r.irate ≈ (8/2)*exp(-Float64(f32(0.5))*13.605692/(0.861707*Tc))  rtol=1e-6
        @test Radix.rate(mk(2, 1), Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels) == r
        # the fit temperature is floored where eij/kT would exceed 50
        cold = Radix.rate(c, Radix.Cell(0.01, 1e3, 2e3, 1e4); levels=levels)
        @test isfinite(cold.irate) && cold.irate > 0
        # records with other than 5 or 9 knots are skipped
        odd = Radix.ElectronCollision(Int8(3), "x", Int32(2), Radix.Transition(Int32(1), Int32(2)),
            Int32(1), Int32(7), f32(0.5), f32(C0), f32.(kn[1:4]))
        @test Radix.rate(odd, Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels).frate == 0
        # a spline that dips below 0 gives zero rates, not negative ones
        neg = Radix.ElectronCollision(Int8(3), "x", Int32(1), Radix.Transition(Int32(1), Int32(2)),
            Int32(1), Int32(7), f32(20.76), f32(1.3), f32[0.0, 0.002592, 0.01144, 0.02149, 0.0363])
        @test Radix.chianti_upsilon(1, 20.76, 1.3, f32[0.0, 0.002592, 0.01144, 0.02149, 0.0363], 1e3) < 0
        @test Radix.rate(neg, Radix.Cell(0.1, 1e3, 2e3, 1e4); levels=levels).irate == 0
        @test Radix.rate(mk(1, 9), Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels).frate == 0
        @test Radix.rate(mk(1, 2; ΔE=0), Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels).frate == 0
    end

    @testset "CollisionProb" begin
        # the numbers are checked against ucalc in ucalc_tests.jl
        none_lv = Radix.level_table(Radix.AtomicLevel[])
        cpx = Radix.CollisionProb(Int8(8), "x", Radix.Transition(Int32(1), Int32(2)), Int32(8), Int32(7))
        @test Radix.rate(cpx, Radix.Cell(1.0, 0.0, 1e4, 1e4); levels=none_lv).frate == 0
    end

    @testset "functor convenience" begin
        c = Radix.RadRecomb(Int32(1), "rr", Int32[1], f32[1e-12, 0.7])
        @test c(cell) == Radix.rate(c, cell)
        @test c(cell; index=true) == Radix.rate(c, cell; index=true)
        cells = [Radix.Cell(T, nh, ne, 1e4) for T in (0.5, 1.0, 2.0)]
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

            # levels: the reals must be read from the right place (FITSFiles
            # work-around in load); every level has a positive weight
            levels = level_table(db)
            @test length(levels) == 38235
            @test all(l -> l.g > 0, values(levels))
            s2 = levels[(122, 1)]                       # S II ground level 3p3 4S
            @test (s2.E, s2.g) == (0, 4)
            @test s2.E_inf ≈ 23.4

            # AtomicLine2: the decay runs downward in energy and all results are
            # finite; the 339 records without a wavelength give nothing, as in ucalc
            cell0 = Radix.Cell(1.0, 1e4, 1e4, 1e4)
            nbad = 0; nord = 0; nlines = 0; nskip = 0
            for r in db
                r isa Radix.AtomicLine2 || continue
                nlines += 1
                x = Radix.rate(r, cell0; levels=levels, mass=16.0)
                x.init == 0 && (r.λ == 0 ? (nskip += 1) : (nbad += 1); continue)
                levels[(r.ion, x.init)].E >= levels[(r.ion, x.final)].E || (nord += 1)
                all(isfinite, (x.frate, x.fenergy, x.opacity)) || (nbad += 1)
            end
            @test nlines == 730369
            @test nskip == 339
            @test nbad == 0
            @test nord == 0

            # ... and with a spectrum: photoexcitation rates are finite and non-negative
            E = [0.1*exp(0.0015*(i - 1)) for i in 1:9999]
            rad = Radiation(E, 1e16 .* E .^ -1.0)
            nbad = 0
            for r in db
                r isa Radix.AtomicLine2 || continue
                x = Radix.rate(r, cell0; levels=levels, mass=16.0, radiation=rad, cfrac=0.25)
                x.init == 0 && continue
                (isfinite(x.irate) && x.irate >= 0 && isfinite(x.ienergy)) || (nbad += 1)
            end
            @test nbad == 0

            # ElectronImpact1: levels always resolve, the 188 degenerate pairs
            # (equal energies) are skipped as in ucalc, and the rates are finite,
            # non-negative and obey detailed balance at 1e3, 1e4 and 1e5 K
            nrec = 0; nnone = 0; nbad = 0
            for T in (0.1, 1.0, 10.0)
                c = Radix.Cell(T, 1e4, 1e4, 1e4)
                for r in db
                    r isa Radix.ElectronImpact1 || continue
                    T == 0.1 && (nrec += 1)
                    x = Radix.rate(r, c; levels=levels)
                    x.init == 0 && (T == 0.1 && (nnone += 1); continue)
                    ok = isfinite(x.frate) && isfinite(x.irate) && x.frate >= 0 && x.irate >= 0
                    lo, up = levels[(r.ion, x.init)], levels[(r.ion, x.final)]
                    ok &= x.irate == 0 ||
                        isapprox(x.frate/x.irate, (up.g/lo.g)*Radix.expo(-(up.E - lo.E)/(0.861707*T)); rtol=1e-6)
                    ok || (nbad += 1)
                end
            end
            @test nrec == 87231
            @test nnone == 188
            @test nbad == 0

            # ElectronCollision: finite and obeying detailed balance with the
            # record's own transition energy. Υ is clamped at 0 (ucalc does not,
            # and its spline dips slightly below 0 for 4 records).
            nrec = 0; nnone = 0; nbad = 0; nneg = 0
            for T in (0.1, 1.0, 10.0)
                c = Radix.Cell(T, 1e4, 1e4, 1e4)
                for r in db
                    r isa Radix.ElectronCollision || continue
                    T == 0.1 && (nrec += 1)
                    x = Radix.rate(r, c; levels=levels)
                    x.init == 0 && (T == 0.1 && (nnone += 1); continue)
                    lo, up = levels[(r.ion, x.init)], levels[(r.ion, x.final)]
                    ok = isfinite(x.frate) && isfinite(x.irate)
                    ok &= x.irate == 0 || isapprox(x.frate/x.irate,
                        (up.g/lo.g)*Radix.expo(-Float64(r.ΔE)*13.605692/(0.861707*T)); rtol=1e-6)
                    ok || (nbad += 1)
                    (x.frate < 0 || x.irate < 0) && (nneg += 1)
                end
            end
            @test nrec == 23232
            @test nnone == 0
            @test nbad == 0
            @test nneg == 0

            # ElectronImpact2 (CHIANTI 2016): finite, non-negative, detailed balance
            nrec = 0; nbad = 0
            for T in (0.1, 1.0, 10.0)
                c = Radix.Cell(T, 1e4, 1e4, 1e4)
                for r in db
                    r isa Radix.ElectronImpact2 || continue
                    T == 0.1 && (nrec += 1)
                    x = Radix.rate(r, c; levels=levels)
                    x.init == 0 && (nbad += 1; continue)
                    lo, up = levels[(r.ion, x.init)], levels[(r.ion, x.final)]
                    ok = isfinite(x.frate) && isfinite(x.irate) && x.frate >= 0 && x.irate >= 0
                    ok &= x.irate == 0 || isapprox(x.frate/x.irate,
                        (up.g/lo.g)*Radix.expo(-Float64(r.ΔE)*13.605692/(0.861707*T)); rtol=1e-6)
                    ok || (nbad += 1)
                end
            end
            @test nrec == 841
            @test nbad == 0

            # CollisionProb: all rates finite and non-negative over four (T, nₑ)
            nbad = 0; nrec = 0
            for (T, ne) in ((1.0, 1e4), (10.0, 1e8), (100.0, 1e10), (1000.0, 1e12))
                c = Radix.Cell(T, 0.0, ne, ne)
                for r in db
                    r isa Radix.CollisionProb || continue
                    T == 1.0 && (nrec += 1)
                    x = Radix.rate(r, c; levels=levels)
                    (isfinite(x.frate) && isfinite(x.irate) && x.frate >= 0 && x.irate >= 0) || (nbad += 1)
                end
            end
            @test nrec == 6015
            @test nbad == 0

            # ParPhotoIonize1/2 against a power-law spectrum: all results finite and
            # photoionization rates non-negative, except for the 18 type-53 records whose
            # tables start hundreds of eV above the threshold (negative in XSTAR too)
            counts = level_counts(levels)
            E = [0.1*exp(0.0015*(i - 1)) for i in 1:9999]
            rad = Radiation(E, 1e16 .* E .^ -1.0)
            nbad = 0; nneg = 0; nev = 0
            for (T, ne) in ((0.3, 5e3), (10.0, 1e10))
                c = Radix.Cell(T, 1e3, ne, 1e10)
                for r in db
                    (r isa Radix.ParPhotoIonize1 || r isa Radix.ParPhotoIonize2) || continue
                    x = Radix.rate(r, c; levels=levels, radiation=rad, nlev=counts[r.ion], lfast=3)
                    x.init == 0 && continue
                    nev += 1
                    all(isfinite, (x.frate, x.irate, x.fenergy, x.ienergy, x.opacity)) || (nbad += 1)
                    (x.frate < 0 || x.irate < 0) && (nneg += 1)
                end
            end
            @test nbad == 0
            @test nneg == 36
            @test nev > 598000

            # the other photoionization types: finite results over the same spectrum
            kinds = Dict(Radix.ParPhotoIonize3 => 0, Radix.PhotoionizeDelta => 0,
                Radix.PhotoionizeDamp => 0, Radix.PhotoionizeFeKedge => 0,
                Radix.PhotoionizeSuper => 0, Radix.PhotoRecombX => 0)
            nbad = 0
            for (T, ne) in ((0.3, 5e3), (10.0, 1e10))
                c = Radix.Cell(T, 1e3, ne, 1e10)
                for r in db
                    typeof(r).name.wrapper in keys(kinds) || continue
                    x = Radix.rate(r, c; levels=levels, radiation=rad, nlev=counts[r.ion])
                    x.init == 0 && continue
                    kinds[typeof(r).name.wrapper] += 1
                    all(isfinite, (x.frate, x.irate)) || (nbad += 1)
                end
            end
            @test nbad == 0
            @test all(>(0), values(kinds))
            @test kinds[Radix.PhotoionizeSuper] > 0 && kinds[Radix.PhotoRecombX] > 0
        end
    else
        @info "Skipping atdb.fits parsing test (set RADIX_ATDB to enable)"
    end
end

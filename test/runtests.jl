using Radix
using Test

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
        # Hand-built level table; the levels are listed upper-first on purpose
        lv(level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1,2,0,1,level,7],
            f32[E, g, 1, 13.6])
        levels = Radix.level_table([lv(1, 0.0, 2), lv(2, 10.2, 8)])
        # ucalc label 50: flin = 1e-16 aij ggup elin²/(0.667274 gglo);
        # vtherm = max(vturb 1e5, 1.3e6/sqrt(a/t)); sigma = 0.02655 flin elin 1e-8/vtherm;
        # ans1 = aij; ans4 = ans1 (12398.54/elin) ergsev
        λ, A = 1215.67, 6.265e8
        c = Radix.AtomicLine2(Int32(4), "lya", Int32[2,1,1,7], f32[λ, 0.4162, A])
        r = Radix.rate(c, cell; levels=levels, mass=1.01, vturb=1.0)
        A64, λ64 = Float64(f32(A)), Float64(f32(λ))      # the data are single precision
        flin = 1e-16*A64*8*λ64^2/(0.667274*2)
        vth = max(1e5, 1.3e6/sqrt(1.01/T))
        @test (r.init, r.final) == (2, 1)          # upper first, whatever the stored order
        @test r.frate == A64
        @test r.irate == 0
        @test r.opacity ≈ 0.02655*flin*λ64*1e-8/vth
        @test r.fenergy ≈ A64*12398.54/λ64*1.602197e-12
        # stored lower-first gives the same answer
        c2 = Radix.AtomicLine2(Int32(4), "lya", Int32[1,2,1,7], f32[λ, 0.4162, A])
        @test Radix.rate(c2, cell; levels=levels, mass=1.01) == Radix.rate(c, cell; levels=levels, mass=1.01)
        # no wavelength or two-photon: no opacity; unknown level: nothing
        c3 = Radix.AtomicLine2(Int32(4), "x", Int32[1,2,1,7], f32[1e9, 0, A])
        @test Radix.rate(c3, cell; levels=levels, mass=1.01).opacity == 0
        c4 = Radix.AtomicLine2(Int32(9), "2ph", Int32[1,2,1,7], f32[λ, 0, A])
        @test Radix.rate(c4, cell; levels=levels, mass=1.01).opacity == 0
        # no wavelength: finite power from the level energy difference
        c6 = Radix.AtomicLine2(Int32(4), "x", Int32[1,2,1,7], f32[0, 0, A])
        @test Radix.rate(c6, cell; levels=levels, mass=1.01).fenergy ≈ A64*Float64(f32(10.2))*1.602197e-12
        c5 = Radix.AtomicLine2(Int32(4), "x", Int32[1,9,1,7], f32[λ, 0, A])
        @test Radix.rate(c5, cell; levels=levels, mass=1.01).frate == 0
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
        # ucalc label 51: tk = max(1e4 t, 2.8777e6/elin), eij = ΔE Ry 13.598,
        #   cji = 8.626e-8 Υ/(sqrt(t) ggup), cij = cji ggup exp(-eij/(0.861707 t))/gglo
        Tc = 20.0
        c = mk(1, 2)
        r = Radix.rate(c, Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels)
        eij = Float64(f32(0.5))*13.598
        tk = max(Tc*1e4, 2.8777e6/(12398.54/eij))
        Υ = Radix.chianti_upsilon(2, Float64(f32(0.5)), Float64(f32(C0)), Float64.(f32.(kn)), tk)
        cji = 8.626e-8*Υ/sqrt(Tc)/8
        cij = cji*8*exp(-eij/(0.861707*Tc))/2
        @test (r.init, r.final) == (1, 2)
        @test r.irate ≈ cji*2e3  rtol=1e-6
        @test r.frate ≈ cij*2e3  rtol=1e-6
        @test Radix.rate(mk(2, 1), Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels) == r
        # the fit temperature is floored where eij/kT would exceed 50
        cold = Radix.rate(c, Radix.Cell(0.01, 1e3, 2e3, 1e4); levels=levels)
        @test isfinite(cold.irate) && cold.irate > 0
        # a spline that dips below 0 gives zero rates, not negative ones
        neg = Radix.ElectronCollision(Int8(3), "x", Int32(1), Radix.Transition(Int32(1), Int32(2)),
            Int32(1), Int32(7), f32(20.76), f32(1.3), f32[0.0, 0.002592, 0.01144, 0.02149, 0.0363])
        @test Radix.chianti_upsilon(1, 20.76, 1.3, f32[0.0, 0.002592, 0.01144, 0.02149, 0.0363], 1e3) < 0
        @test Radix.rate(neg, Radix.Cell(0.1, 1e3, 2e3, 1e4); levels=levels).irate == 0
        @test Radix.rate(mk(1, 9), Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels).frate == 0
        @test Radix.rate(mk(1, 2; ΔE=0), Radix.Cell(Tc, 1e3, 2e3, 1e4); levels=levels).frate == 0
    end

    @testset "CollisionProb" begin
        # against XSTAR's own Fortran (ucalc 63, double-precision build), 800 cases
        # from real records; the file lists ni li nf lf Z T E1 E2 g1 g2 ne ans1 ans2
        nz = 0
        for line in eachline(joinpath(@__DIR__, "reference", "collisionprob_cases.txt"))
            startswith(line, "#") && continue
            ni, li, nf, lf, Zp, Tk, e1, e2, g1, g2, nel, a1, a2 = parse.(Float64, split(line))
            mk(level, n, L, E, g) = Radix.AtomicLevel(Int32(13), "",
                Int32[n, 2, L, Zp, level, 7], f32[E, g, 1, 13.6])
            lv = Radix.level_table([mk(1, ni, li, e1, g1), mk(2, nf, lf, e2, g2)])
            cp = Radix.CollisionProb(Int8(8), "x", Radix.Transition(Int32(1), Int32(2)),
                Int32(Zp), Int32(7))
            r = Radix.rate(cp, Radix.Cell(Tk/1e4, 0.0, nel, nel); levels=lv)
            @test r.frate ≈ a1  rtol=1e-6 atol=1e-300
            @test r.irate ≈ a2  rtol=1e-6 atol=1e-300
            nz += (a1 != 0)
        end
        @test nz > 100                            # the reference is not mostly zeros
        # unknown level: nothing
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

            # AtomicLine2: every record resolves to two known levels, the decay
            # runs downward in energy, and all results are finite
            cell0 = Radix.Cell(1.0, 1e4, 1e4, 1e4)
            nbad = 0; nord = 0; nlines = 0
            for r in db
                r isa Radix.AtomicLine2 || continue
                nlines += 1
                x = Radix.rate(r, cell0; levels=levels, mass=16.0)
                x.init == 0 && (nbad += 1; continue)
                levels[(r.ion, x.init)].E >= levels[(r.ion, x.final)].E || (nord += 1)
                all(isfinite, (x.frate, x.fenergy, x.opacity)) || (nbad += 1)
            end
            @test nlines == 730369
            @test nbad == 0
            @test nord == 0

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
                        (up.g/lo.g)*Radix.expo(-Float64(r.ΔE)*13.598/(0.861707*T)); rtol=1e-6)
                    ok || (nbad += 1)
                    (x.frate < 0 || x.irate < 0) && (nneg += 1)
                end
            end
            @test nrec == 23232
            @test nnone == 0
            @test nbad == 0
            @test nneg == 0

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
        end
    else
        @info "Skipping atdb.fits parsing test (set RADIX_ATDB to enable)"
    end
end

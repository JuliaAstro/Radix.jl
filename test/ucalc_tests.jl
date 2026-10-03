# Radix rates against XSTAR's real ucalc (Fortran 90, ftools/xstar/xstarlib/src).
# test/reference/ucalc/data/typeNN.in holds cases sampled from atdb.fits (records,
# conditions and the level data they need) and typeNN.out what ucalc returns for
# them (see test/reference/ucalc/README.md). The driver is built like HEASoft's
# xstar, so XSTAR's single-precision literals limit the agreement to ~1e-6.

const UCALC_DATA = joinpath(@__DIR__, "reference", "ucalc", "data")

# outputs of the driver: ans1..ans6, idest1..idest4, opakab
function ucalc_cases(name)
    L = readlines(joinpath(UCALC_DATA, name * ".in"))
    O = [parse.(Float64, split(l)) for l in readlines(joinpath(UCALC_DATA, name * ".out"))]
    cases = NamedTuple[]
    i = 1
    while i <= length(L)
        ndesc, rtype, nrdt, nidt = parse.(Int, split(L[i]))
        ints = parse.(Int, split(L[i+1]))
        reals = parse.(Float64, split(L[i+2]))
        cond = parse.(Float64, split(L[i+3]))
        nlev, nspec = parse.(Int, split(L[i+4]))
        lv = [parse.(Float64, split(L[i+4+k])) for k in 1:nspec]
        push!(cases, (; ndesc, rtype, ints, reals, cond, nlev, lv))
        i += 5 + nspec
    end
    @assert length(cases) == length(O)
    zip(cases, O)
end

function ucalc_inputs(c)
    coef = Radix.ratemap[c.ndesc](Int32(c.rtype), "", Int32.(c.ints), Float32.(c.reals))
    t, xpx, xee, xh0 = c.cond[1:4]
    cell = Radix.Cell(t, xh0, xpx*xee, xpx)
    levels = Radix.level_table([Radix.AtomicLevel(Int32(13), "",
        Int32[v[5], v[6], v[7], 1, v[1], c.ints[end]], Float32[v[2], v[3], 1, v[4]])
        for v in c.lv])
    (coef, cell, levels)
end

# `expected(o)` maps the ucalc outputs onto the named fields Radix returns
function check_ucalc(name; call, expected, rtol=1e-5)
    n = 0
    for (c, o) in ucalc_cases(name)
        coef, cell, levels = ucalc_inputs(c)
        got = call(coef, cell, levels, c)
        for (k, v) in pairs(expected(o))
            # a skipped rate has no transition in Radix, ucalc still reports its levels
            k in (:init, :final) && got.init == 0 && continue
            @test isapprox(getfield(got, k), v; rtol=rtol, atol=0)
        end
        n += 1
    end
    n
end

level_call(coef, cell, levels, c) = Radix.rate(coef, cell; levels=levels)
plain_call(coef, cell, levels, c) = Radix.rate(coef, cell)
direct(o; swap=false) = (; frate=o[1], irate=o[2], init=swap ? o[8] : o[7], final=swap ? o[7] : o[8])

@testset "Radix rates vs XSTAR ucalc" begin
    @testset "type 1 RadRecomb" begin
        @test check_ucalc("type01"; call=plain_call, expected=direct) == 100
    end
    @testset "type 2 ChargeExH0" begin
        @test check_ucalc("type02";
            call=(co, ce, lv, c) -> Radix.rate(co, ce; nlev=c.nlev), expected=direct) == 100
    end
    @testset "type 9 ChargeExHe" begin
        @test check_ucalc("type09";
            call=(co, ce, lv, c) -> Radix.rate(co, ce; nlev=c.nlev), expected=direct) == 100
    end
    @testset "type 30 TotRadRecombH" begin
        @test check_ucalc("type30"; call=plain_call, expected=direct) == 100
    end
    @testset "type 38 TotRadRecomb" begin
        @test check_ucalc("type38"; call=plain_call, expected=direct) == 100
    end
    @testset "type 50 AtomicLine2 (without photoexcitation)" begin
        # ucalc returns the photoexcitation in ans1 (0 here: the driver's spectrum is
        # empty), the decay in ans2, minus the emitted power in ans3, the opacity last
        @test check_ucalc("type50";
            call=(co, ce, lv, c) -> Radix.rate(co, ce; levels=lv, mass=c.cond[13],
                vturb=c.cond[6], pesc=c.cond[11] + c.cond[12], nlev=c.nlev),
            expected=o -> (; frate=o[2], irate=o[1], fenergy=-o[3], ienergy=-o[4],
                init=o[7], final=o[8], opacity=o[11]), rtol=1e-6) == 100
    end
    # ElectronCollision returns the lower level first and clamps Υ at 0; ucalc does not
    clamp0(o) = (; frate=max(0, o[1]), irate=max(0, o[2]), init=o[8], final=o[7])
    @testset "type 51 ElectronCollision" begin
        @test check_ucalc("type51"; call=level_call, expected=clamp0) == 100
        @test check_ucalc("type51n9"; call=level_call, expected=clamp0) == 16   # 9-point fits
    end
    @testset "type 56 ElectronImpact1" begin
        @test check_ucalc("type56"; call=level_call, expected=direct) == 100
    end
    @testset "type 63 CollisionProb" begin
        @test check_ucalc("type63"; call=level_call, expected=direct) == 200
    end
    @testset "type 98 ElectronImpact2" begin
        @test check_ucalc("type98"; call=level_call, expected=direct) == 100
    end
end

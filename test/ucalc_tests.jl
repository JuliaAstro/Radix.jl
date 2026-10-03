# Radix rates against XSTAR's real ucalc (Fortran 90, ftools/xstar/xstarlib/src).
# test/reference/ucalc/data/typeNN.in holds cases sampled from atdb.fits (records,
# conditions and the level data they need) and typeNN.out what ucalc returns for
# them (see test/reference/ucalc/README.md). The driver is built like HEASoft's
# xstar, so XSTAR's single-precision literals limit the agreement to ~1e-6.

const UCALC_DATA = joinpath(@__DIR__, "reference", "ucalc", "data")

# Fortran writes exponents of three digits without the E (1.5-133)
fortran_float(s) = parse(Float64, replace(s, r"(\d)([+-]\d{3})$" => s"\1E\2"))

# outputs of the driver: ans1..ans6, idest1..idest4, opakab, four checksums
function ucalc_cases(name)
    L = readlines(joinpath(UCALC_DATA, name * ".in"))
    O = [fortran_float.(split(l)) for l in readlines(joinpath(UCALC_DATA, name * ".out"))]
    cases = NamedTuple[]
    i = 1
    while i <= length(L)
        ndesc, rtype, nrdt, nidt = parse.(Int, split(L[i]))
        ints = parse.(Int, split(L[i+1]))
        reals = parse.(Float64, split(L[i+2]))
        cond = parse.(Float64, split(L[i+3]))
        rad = parse.(Float64, split(L[i+4]))               # e0 dlnE index norm kpar gpar epar lfast
        nlev, nspec = parse.(Int, split(L[i+5]))
        lv = [parse.(Float64, split(L[i+5+k])) for k in 1:nspec]
        push!(cases, (; ndesc, rtype, ints, reals, cond, rad, nlev, lv))
        i += 6 + nspec
    end
    @assert length(cases) == length(O)
    zip(cases, O)
end

# the driver's radiation field: epi(i) = e0 exp(dlnE (i-1)), bremsa = norm epi^-index, 9999 points
function ucalc_radiation(c)
    e0, dlnE, index, norm = c.rad[1:4]
    E = [e0*exp(dlnE*(i - 1)) for i in 1:9999]
    Radix.Radiation(E, norm .* E .^ (-index))
end

# the parent ion of the photoionization records (position among the integers)
parent_ion(c) = c.ints[c.ndesc == 59 ? 5 : c.ndesc in (70, 99) ? 9 : 6]
parent_level(c) = c.ints[c.ndesc == 59 ? 4 : c.ndesc in (70, 99) ? 8 : 5]

function ucalc_inputs(c)
    coef = Radix.ratemap[c.ndesc](Int32(c.rtype), "", Int32.(c.ints), Float32.(c.reals))
    t, xpx, xee, xh0 = c.cond[1:4]
    cell = Radix.Cell(t, xh0, xpx*xee, xpx)
    recs = [Radix.AtomicLevel(Int32(13), "",
        Int32[v[5], v[6], v[7], 1, v[1], c.ints[end]], Float32[v[2], v[3], 1, v[4]])
        for v in c.lv]
    # photoionization leaving an excited level of the next ion: its weight and energy (the
    # driver gets them as kpar, gpar, epar)
    if c.ndesc in (49, 53, 59, 70, 99) && c.rad[5] > 1
        push!(recs, Radix.AtomicLevel(Int32(13), "", Int32[1, 2, 0, 1, c.rad[5], parent_ion(c)],
            Float32[c.rad[7], c.rad[6], 1, 0]))
    end
    (coef, cell, Radix.level_table(recs))
end

# `expected(o)` maps the ucalc outputs onto the named fields Radix returns
function check_ucalc(name; call, expected, rtol=1e-5, atol=0)
    n = 0
    for (c, o) in ucalc_cases(name)
        coef, cell, levels = ucalc_inputs(c)
        got = call(coef, cell, levels, c)
        for (k, v) in pairs(expected(o, c))
            # a skipped rate has no transition in Radix, ucalc still reports its levels
            k in (:init, :final) && got.init == 0 && continue
            @test isapprox(getfield(got, k), v; rtol=rtol, atol=atol)
        end
        n += 1
    end
    n
end

level_call(coef, cell, levels, c) = Radix.rate(coef, cell; levels=levels)
plain_call(coef, cell, levels, c) = Radix.rate(coef, cell)
direct(o, c) = (; frate=o[1], irate=o[2], init=o[7], final=o[8])
# the collision types also return the energies rate × ΔE (ans6 forward, ans5 inverse)
with_energy(o, c) = (; direct(o, c)..., fenergy=o[6], ienergy=o[5])

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
    @testset "type 50 AtomicLine2" begin
        # ucalc returns the photoexcitation in ans1, the decay in ans2, minus the energy of
        # the decays in ans3 and of the photoexcitations in ans4, and the opacity last
        @test check_ucalc("type50";
            call=(co, ce, lv, c) -> Radix.rate(co, ce; levels=lv, mass=c.cond[13],
                vturb=c.cond[6], pesc=c.cond[11] + c.cond[12], nlev=c.nlev,
                radiation=ucalc_radiation(c), cfrac=c.cond[8]),
            expected=(o, c) -> (; frate=o[2], irate=o[1], fenergy=-o[3], ienergy=-o[4],
                init=o[7], final=o[8], opacity=o[11]), rtol=1e-6) == 100
    end
    # ElectronCollision returns the lower level first and clamps Υ at 0; ucalc does not
    clamp0(o, c) = (; frate=max(0, o[1]), irate=max(0, o[2]), init=o[8], final=o[7],
        fenergy=max(0, o[6]), ienergy=max(0, o[5]))
    @testset "type 51 ElectronCollision" begin
        @test check_ucalc("type51"; call=level_call, expected=clamp0) == 100
        @test check_ucalc("type51n9"; call=level_call, expected=clamp0) == 16   # 9-point fits
    end
    @testset "type 56 ElectronImpact1" begin
        @test check_ucalc("type56"; call=level_call, expected=with_energy) == 100
    end
    @testset "type 63 CollisionProb" begin
        @test check_ucalc("type63"; call=level_call, expected=with_energy) == 200
    end
    # photoionization (values below 1e-100 are physically zero and depend on the last digit of
    # single-precision constants in the exponent): also the opacity and recombination-emissivity arrays the integrals add
    # to; the energies measured from the levels (fenergy2, ienergy2) are only compared for
    # parent level 1, since ucalc reads stale memory for the energy of higher ones
    photo_call = (co, ce, lv, c; kw...) -> begin
        op = Radix.Opacity(9999)
        r = Radix.rate(co, ce; levels=lv, radiation=ucalc_radiation(c), nlev=c.nlev,
            ptmp=(c.cond[11], c.cond[12]), abund=(c.cond[9], c.cond[10]),
            lfast=Int(c.rad[8]), opacity=op, kw...)
        merge(r, (; sum_total=sum(op.total), sum_continuum=sum(op.continuum),
            sum_em1=sum(op.emissivity[1, :]), sum_em2=sum(op.emissivity[2, :])))
    end
    photo_expected = (o, c) -> merge((; frate=o[1], irate=o[2], ienergy=-o[3], fenergy=-o[4],
            init=o[7], final=o[8], opacity=o[11], sum_total=o[12], sum_continuum=o[13],
            sum_em1=o[14], sum_em2=o[15]),
        parent_level(c) == 1 || c.ndesc == 59 ? (; ienergy2=-o[5], fenergy2=-o[6]) : (;))
    @testset "type 49 ParPhotoIonize1" begin
        @test check_ucalc("type49"; call=photo_call, expected=photo_expected, rtol=5e-6, atol=1e-100) == 96
    end
    @testset "type 53 ParPhotoIonize2" begin
        @test check_ucalc("type53"; call=photo_call, expected=photo_expected, rtol=5e-6, atol=1e-100) == 96
    end
    @testset "type 59 ParPhotoIonize3" begin
        @test check_ucalc("type59"; call=photo_call, expected=photo_expected, rtol=1e-6, atol=1e-100) == 100
    end
    @testset "type 70 PhotoionizeSuper" begin
        # ucalc's jkion = 1: the density is limited to 1e8
        @test check_ucalc("type70"; call=(co, ce, lv, c) -> photo_call(co, ce, lv, c; neutral=true),
            expected=photo_expected, rtol=1e-6, atol=1e-100) == 80
    end
    @testset "type 99 PhotoRecombX" begin
        # the energy correction is a difference of nearly equal terms (1.1e-6 at worst)
        @test check_ucalc("type99"; call=photo_call, expected=photo_expected, rtol=5e-6, atol=1e-100) == 120
        # an excited level of the parent ion (made up): the energies from the levels are stale in ucalc
        @test check_ucalc("type99par"; call=photo_call, expected=photo_expected, rtol=5e-6, atol=1e-100) == 40
    end
    @testset "type 74 PhotoionizeDelta" begin
        @test check_ucalc("type74";
            call=(co, ce, lv, c) -> Radix.rate(co, ce; levels=lv, radiation=ucalc_radiation(c),
                nlev=c.nlev),
            expected=(o, c) -> (; frate=o[1], irate=o[2], init=o[7], final=o[8]),
            rtol=5e-6, atol=1e-100) == 100
    end
    @testset "type 85 PhotoionizeFeKedge" begin
        # recombination and the opacity are zeroed in ucalc; the opacity arrays are filled
        @test check_ucalc("type85"; call=photo_call, expected=photo_expected, rtol=1e-6, atol=1e-100) == 100
    end
    @testset "type 88 PhotoionizeDamp" begin
        # only the photoionization rate and the opacity are returned; the arrays are still filled.
        # The opacity is a difference of two nearly equal terms, so it magnifies the
        # single-precision constants of XSTAR (2.5e-5 at worst), hence the tolerance
        @test check_ucalc("type88"; call=photo_call, expected=photo_expected, rtol=1e-4, atol=1e-100) == 48
    end
    @testset "type 53 with tables that start far above the threshold" begin
        # these records give negative photoionization rates in XSTAR too
        @test check_ucalc("type53odd"; call=photo_call, expected=photo_expected, rtol=5e-5,
            atol=1e-100) == 12
    end
    @testset "type 54 RadiativeProb" begin
        @test check_ucalc("type54"; call=level_call,
            expected=(o, c) -> (; frate=o[1], irate=o[2], ienergy=-o[3], fenergy=-o[4], init=o[7], final=o[8])) == 240
    end
    @testset "type 57 EffectiveCharge" begin
        # collisional ionization and three-body recombination; ucalc gives rates for 180 of the 1600 sampled cases
        @test check_ucalc("type57";
            call=(co, ce, lv, c) -> Radix.rate(co, ce; levels=lv, nlev=c.nlev),
            expected=(o, c) -> (; frate=o[1], irate=o[2], fenergy=o[6], ienergy=o[5], init=o[7], final=o[8]),
            rtol=2e-6, atol=1e-100) == 240
    end
    @testset "type 77 CollisionSuper" begin
        @test check_ucalc("type77"; call=(co, ce, lv, c) -> Radix.rate(co, ce; levels=lv, nlev=c.nlev),
            expected=with_energy, rtol=1e-6, atol=1e-100) == 240
    end
    @testset "type 86 IronKAuger" begin
        @test check_ucalc("type86"; call=(co, ce, lv, c) -> Radix.rate(co, ce; nlev=c.nlev), expected=direct) == 160
    end
    @testset "type 71 RadiativeSuper" begin
        @test check_ucalc("type71";
            call=(co, ce, lv, c) -> Radix.rate(co, ce; levels=lv, nlev=c.nlev, mass=c.cond[13],
                vturb=c.cond[6], ptmp=(c.cond[11], c.cond[12])),
            expected=(o, c) -> (; frate=o[1], irate=o[2], ienergy=-o[3], fenergy=-o[4], init=o[7], final=o[8],
                opacity=o[11]),
            rtol=1e-6, atol=1e-100) == 240
        # calcium I and II (ions 96 and 97), whose decay rate is capped at 1e10
        @test check_ucalc("type71ions";
            call=(co, ce, lv, c) -> Radix.rate(co, ce; levels=lv, nlev=c.nlev, mass=c.cond[13],
                vturb=c.cond[6], ptmp=(c.cond[11], c.cond[12])),
            expected=(o, c) -> (; irate=o[2], ienergy=-o[3], opacity=o[11]),
            rtol=1e-6, atol=1e-100) == 8
    end
    @testset "type 95 CollisionIonize" begin
        @test check_ucalc("type95"; call=(co, ce, lv, c) -> Radix.rate(co, ce; levels=lv, nlev=c.nlev),
            expected=with_energy, rtol=1e-6, atol=1e-100) == 160
    end
    @testset "type 98 ElectronImpact2" begin
        @test check_ucalc("type98"; call=level_call, expected=with_energy) == 100
    end
end

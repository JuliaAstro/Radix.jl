# The rates are differentiable with ForwardDiff: the derivatives of `frate` and `irate` with
# respect to the temperature, and to the density (scaling nₕ, nₑ and ntot together), agree with
# central finite differences on the cases of the ucalc fixtures (ucalc_tests.jl).
# The conditions are moved off the table nodes, where the interpolations have kinks. The finite
# differences use the same dual-number code path (values only), which for type 99 avoids the
# single-precision path that plain floats take.

using ForwardDiff

const AD_STEP = 1e-4                  # relative step of the finite differences
const AD_RTOL = 2e-3
# Type 70 integrates with a step that adapts to a 1% convergence test and sums the Milne
# relation until it changes by less than 1%, which makes it slightly non-smooth: its finite
# differences are noisy at the 1% level while the derivative stays within it.
const AD_RTOLS = Dict("type70" => 3e-2)
# the fewest derivatives a case must compare (type 22 has three records, half of them above its 6e4 K limit; type 10 has four)
const AD_MIN_CHECKS = Dict("type22" => 4, "type10" => 4)
const AD_OFFSET = 1.2345              # moves T and the densities off the table nodes

ad_photo(co, ce, lv, c; kw...) = Radix.rate(co, ce; radiation=ucalc_radiation(c), ptmp=(c.cond[11], c.cond[12]), abund=(c.cond[9], c.cond[10]), lfast=Int(c.rad[8]), kw...)
ad_with_opacity(co, ce, lv, c) = ad_photo(co, ce, lv, c; opacity=Radix.Opacity(typeof(ce.T), 9999))

const AD_CASES = [
    "type01" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type02" => (co, ce, lv, c) -> Radix.rate(co, ce; nlev=c.nlev),
    "type08" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type07" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type22" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type39" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type10" => (co, ce, lv, c) -> Radix.rate(co, ce; nlev=c.nlev),   # (the H⁺ density follows the cell)
    "type09" => (co, ce, lv, c) -> Radix.rate(co, ce; nlev=c.nlev),
    "type30" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type38" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type50" => (co, ce, lv, c) -> Radix.rate(co, ce; vturb=c.cond[6],
        pesc=c.cond[11] + c.cond[12], radiation=ucalc_radiation(c), cfrac=c.cond[8]),
    "type91" => (co, ce, lv, c) -> Radix.rate(co, ce; vturb=c.cond[6],
        pesc=c.cond[11] + c.cond[12], radiation=ucalc_radiation(c), cfrac=c.cond[8]),
    "type51" => level_call,
    "type54" => level_call,
    "type56" => level_call,
    "type57" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type60" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type62" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type63" => level_call,
    "type71" => (co, ce, lv, c) -> Radix.rate(co, ce; vturb=c.cond[6], ptmp=(c.cond[11], c.cond[12])),
    "type77" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type76" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type72" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type75" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type73" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type81" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type82" => (co, ce, lv, c) -> Radix.rate(co, ce; vturb=c.cond[6],
        pesc=c.cond[11] + c.cond[12], radiation=ucalc_radiation(c)),
    "type86" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type92hot" => level_call,
    "type95" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type98" => level_call,
    "type49" => ad_photo, "type53" => ad_photo, "type59" => ad_photo,
    "type74" => (co, ce, lv, c) -> Radix.rate(co, ce; radiation=ucalc_radiation(c)),
    "type85" => ad_photo, "type88" => ad_photo, "type99" => ad_photo,
    "type66" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type68" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type69" => (co, ce, lv, c) -> Radix.rate(co, ce),
    "type70" => (co, ce, lv, c) -> ad_photo(co, ce, lv, c; neutral=true),
    "type49 (opacity arrays)" => ad_with_opacity,
]

# the cell with the variable x: `:T` is the temperature, `:n` scales the densities
function ad_cell(ce, variable, x)
    R = typeof(x)
    variable == :T ? Radix.Cell{R}(x, R(ce.nₕ), R(ce.nₑ), R(ce.ntot)) :
        Radix.Cell{R}(R(ce.T), x*ce.nₕ, x*ce.nₑ, x*ce.ntot)
end

# returns the number of derivatives compared and the number that disagree
function ad_check(name, call, variable)
    fixture = first(split(name))
    rtol = get(AD_RTOLS, fixture, AD_RTOL)
    nchecked = nbad = 0
    for (c, _) in ucalc_cases(fixture)
        coef, cell, lv = ucalc_inputs(c)
        x0 = variable == :T ? cell.T*AD_OFFSET : AD_OFFSET
        for k in (:frate, :irate)
            f(x) = getfield(call(coef, ad_cell(cell, variable, x), lv, c), k)
            v = f(x0)
            (v == 0 || !isfinite(v)) && continue
            d = ForwardDiff.derivative(f, x0)
            value(x) = ForwardDiff.value(f(ForwardDiff.Dual{Nothing}(x, 0.0)))
            h = AD_STEP*x0
            fd = (value(x0 + h) - value(x0 - h))/(2h)
            # a derivative that is zero to rounding is compared with the size of the rate
            scale = max(abs(d), abs(fd))
            nchecked += 1
            scale < 1e-8*abs(v)*abs(x0)^-1 || abs(d - fd) <= rtol*scale || (nbad += 1)
        end
    end
    (nchecked, nbad)
end

@testset "ForwardDiff through the rates" begin
    for (name, call) in AD_CASES, variable in (:T, :n)
        @testset "$name, d/d$variable" begin
            nchecked, nbad = ad_check(name, call, variable)
            @test nchecked >= get(AD_MIN_CHECKS, first(split(name)), 21)
            @test nbad == 0
        end
    end
end

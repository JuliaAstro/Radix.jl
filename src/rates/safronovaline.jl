# Radiative transition probabilities of Safronova (Kato et al. 1997): r1 = λ (Å); r2 = not used by ucalc; r3 = A (s⁻¹); i1= i (lower level); i2= k (upper level); i3= Z; i4= ionN

# XSTAR data type: 89

const SafronovaLineDesc = "saf line rad. rates"

struct SafronovaLine{I, R, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # the two levels, in stored order
    Z::I      # atomic number
    ion::I    # ion index (XSTAR ionN)
    λ::R      # wavelength (Å)
    gf::R     # not used
    A::R      # Einstein A (s⁻¹)
    levels::L              # the level data of the database (a Levels)
end

function SafronovaLine(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    SafronovaLine(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., rvec[1:3]..., levels)
end

const safronova_min_wavelength = 1e-48          # Å; a record below it has the decay rate and nothing else
const safronova_f_coeff = Float64(0.667274f0)   # f = 10⁻¹⁶ A g_u λ²/(0.667274 g_l)
const safronova_no_wavelength = Float64(0.99f9) # Å; the lines beyond have no opacity

"""
    rate(coef::SafronovaLine, cell; mass=weight(coef.levels, coef.ion), vturb=1.0, ptmp=(0.5, 0.5), index=false)

Radiative decay of a line from the data of Safronova (XSTAR ucalc type 89): `frate` is the decay rate `A (ptmp₁ + ptmp₂)` from the upper to the lower level (`init` is the upper level by energy; levels `1` to `nlev - 1`)
and `irate` is 0 (ucalc turned the photoexcitation off). With a wavelength above `10⁻⁴⁸` Å `fenergy` is minus the rate times the energy `hc/λ` (ucalc's `ans4` is the positive one) and `opacity` is the line-centre cross section
`0.02655 f λ/v` with `f = 10⁻¹⁶ A g_u λ²/(0.667274 g_l)` and the thermal and turbulent speed `v` of the mass of the element, `mass`, and `vturb` (km/s), without the abundance of the lower level; it is 0 for a wavelength of `0.99×10⁹` Å or more.
The line goes in the continuum bins as the lines of type 50 do not: ucalc's call of `linopac` (with the first line of the ion for the wavelength) is not reproduced.
"""
function rate(coef::SafronovaLine, cell::Cell; mass=weight(coef.levels, coef.ion), vturb=1.0, ptmp=(0.5, 0.5), index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()
    none = (; init=Int32(0), final=Int32(0), frate=0., irate=0., fenergy=0., ienergy=0., opacity=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 >= nlev || i2 <= 0 || i2 >= nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    up, lo = a.E < b.E ? (b, a) : (a, b)
    index && return (; none..., init=up.level, final=lo.level)
    A = Float64(coef.A)
    frate = A*(ptmp[1] + ptmp[2])
    elin = abs(Float64(coef.λ))
    elin <= safronova_min_wavelength && return (; none..., init=up.level, final=lo.level, frate)
    flin = 1e-16*A*Float64(up.g)*elin*elin/(safronova_f_coeff*Float64(lo.g))
    vtherm = sqrt((vturb*cm_per_km)^2 + (K.thermal_speed/sqrt(mass/cell.T))^2)
    opacity = elin > safronova_no_wavelength ? 0.0 : K.line_xsec*flin*elin*cm_per_Å/vtherm
    ener = K.hc_eVÅ_single/elin
    (; init=up.level, final=lo.level, frate, irate=0., fenergy=-frate*ener*K.ergsev, ienergy=0., opacity)
end

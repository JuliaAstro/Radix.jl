# Fluorescence lines: r1 = λ (Å); r2 = f (oscillator strength); r3, r4 = not used; r5 = atomic mass (amu); i1= i (lower level); i2= k (upper level); i3= Z; i4= ionN

# XSTAR data type: 79

const FluorescenceLineDesc = "fluorescence lines"

struct FluorescenceLine{I, R, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # the two levels, in stored order
    Z::I      # atomic number
    ion::I    # ion index (XSTAR ionN)
    λ::R      # wavelength (Å)
    f::R      # oscillator strength
    mass::R   # atomic mass (amu): the fifth real
    levels::L              # the level data of the database (a Levels)
end

function FluorescenceLine(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    FluorescenceLine(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., rvec[1], rvec[2], rvec[5], levels)
end

const fluor_A_coeff = Float64(6.67f7)       # A = 6.67×10⁷ g_l f/(g_u λ²) with λ in μm
const fluor_f_min = Float64(1.01f-12)       # the lines with a smaller oscillator strength or a wavelength above fluor_λ_max have A = fluor_A_fudge
const fluor_λ_max = Float64(1f9)
const fluor_A_fudge = Float64(1f5)          # "a fudge to avoid bad numerics from fine structure"
const fluor_cm_per_μm = 1e-4

"""
    rate(coef::FluorescenceLine, cell; vturb=1.0, ptmp=(0.5, 0.5), index=false)

Radiative decay of a fluorescence line (XSTAR ucalc type 79): `frate` is the rate `A (ptmp₁ + ptmp₂)` from the upper to the lower level (by energy: `init` is the upper level), from the oscillator
strength `f` and the wavelength `λ` (Å) of the record, `A = 6.67×10⁷ g_l f/(g_u λ²)` (λ in μm), or 10⁵ s⁻¹ if `f ≤ 1.01×10⁻¹²` or `λ ≥ 10⁹ Å`; `irate` is 0. `fenergy` is minus the rate times the energy `hc/λ`
(the `ans4` of ucalc is the positive one) and `ienergy` is 0. `opacity` is the line-centre cross section `0.02655 f λ/v` for the thermal and turbulent speed
`v` of the mass of the record (its fifth real) and `vturb` (km/s), without the abundance of the lower level. Both levels must be in `1:nlev`.
"""
function rate(coef::FluorescenceLine, cell::Cell; vturb=1.0, ptmp=(0.5, 0.5), index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()
    none = (; init=Int32(0), final=Int32(0), frate=0., irate=0., fenergy=0., ienergy=0., opacity=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 > nlev || i2 <= 0 || i2 > nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    up, lo = a.E < b.E ? (b, a) : (a, b)
    index && return (; none..., init=up.level, final=lo.level)
    elin = abs(Float64(coef.λ))
    flin = Float64(coef.f)
    elammu = elin*fluor_cm_per_μm
    aij = fluor_A_coeff*Float64(lo.g)*flin/Float64(up.g)/elammu/elammu
    (flin <= fluor_f_min || elin >= fluor_λ_max) && (aij = fluor_A_fudge)
    frate = aij*(ptmp[1] + ptmp[2])
    vtherm = sqrt((vturb*cm_per_km)^2 + (K.thermal_speed/sqrt(Float64(coef.mass)/cell.T))^2)
    opacity = K.line_xsec*flin*elin*cm_per_Å/vtherm
    (; init=up.level, final=lo.level, frate, irate=0., fenergy=-frate*K.ergsev*K.hc_eVÅ_single/elin, ienergy=0., opacity)
end

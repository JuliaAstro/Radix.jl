# Effective collision strengths of He-like ions from Keenan, McCann and Kingston (1987), eq. (2): a quadratic in log₁₀ T
# r1−r3 = coefficients (the first is also the wavelength of the transition in ucalc); i1= i (lower level); i2= j (upper level); i3= Z; i4= ionN

# XSTAR data type: 67

const CollisionKeenanDesc = "Effective collision strengths from Keenan et al."

struct CollisionKeenan{I, R, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # the two levels (ordered by energy for the rates)
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    coeffs::Vector{R}   # coefficients of Υ = c₁ + c₂ log₁₀ T + c₃ log₁₀² T
    levels::L              # the level data of the database (a Levels)
end

function CollisionKeenan(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionKeenan(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., Vector(rvec), levels)
end

const keenan_wavelength_min = 1e-24      # the records with a wavelength (the first coefficient) below this give no rate

"""
    rate(coef::CollisionKeenan, cell; index=false)

Electron-impact excitation and de-excitation of a He-like ion from the quadratic fit of Keenan, McCann and Kingston to the effective collision strength (XSTAR ucalc type 67, `calt67`):
Υ = c₁ + c₂ log₁₀ T + c₃ log₁₀² T with T the temperature in K, at least `T_floor_coeff/|c₁|`, and limited to 0 from below. `irate` is the de-excitation rate `nₑ C Υ/(√T g_u)` and
`frate = irate g_u e^{-Δ/kT}/g_l` (T in 10⁴ K), the exponent `12398.4016/|c₁|/kT` taking the first coefficient as the wavelength in Å. `fenergy` and `ienergy` are the rates times `|c₁|` (ucalc takes the
wavelength for the energy: `dele = elin`). `init` is the lower and `final` the upper level by energy, both in `1:nlev`; a record with `|c₁| ≤ 10⁻²⁴` gives no rates.
"""
function rate(coef::CollisionKeenan, cell::Cell; index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()
    none = (; init=Int32(0), final=Int32(0), frate=0., irate=0., fenergy=0., ienergy=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 > nlev || i2 <= 0 || i2 > nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    lo, up = b.E < a.E ? (b, a) : (a, b)
    index && return (; none..., init=lo.level, final=up.level)
    c1 = Float64(coef.coeffs[1])
    elin = abs(c1)
    elin <= keenan_wavelength_min && return none
    T = max(cell.T*T_unit, K.T_floor_coeff/elin)
    tp = log10(T)
    Υ = max(c1 + Float64(coef.coeffs[2])*tp + Float64(coef.coeffs[3])*tp*tp, 0.0)
    gu, gl = Float64(up.g), Float64(lo.g)
    cji = K.collision_rate_coeff*Υ/sqrt(cell.T)/gu
    cij = cji*gu*expo(-K.hc_eVÅ_single/elin/(K.kT_eV*cell.T))/gl
    frate, irate = cij*cell.nₑ, cji*cell.nₑ
    (; init=lo.level, final=up.level, frate, irate, fenergy=frate*elin*K.ergsev, ienergy=irate*elin*K.ergsev)
end

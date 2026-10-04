# Decay rates for Fe UTA [24]: r1 = λ(Å); r2 = E(k) (eV); r3= g f (i,k);
# r4 = Ar(k,i) (s−1); r5 = Aa(k,i) (s−1); i1= i (lower level);
# i2= k (upper level); i3= ionN

# XSTAR data type: 82

const RadiativeFeDecayDesc = "Fe UTA rad rates"

@with_levels struct RadiativeFeDecay{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    ion::I     # ion index (XSTAR ionN)
    λ::R       # wavelength (Å)
    E::R       # eV
    gf::R      # weighted oscillator strength
    A_rad::R   # s⁻¹
    A_auto::R  # s⁻¹
end

function RadiativeFeDecay(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    RadiativeFeDecay(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3], rvec...)
end

"""
    rate(coef::RadiativeFeDecay, cell; vturb=1.0, pesc=1.0, radiation=nothing, index=false)

Radiative decay and photoexcitation of an Fe UTA line (XSTAR ucalc type 82), like `rate(::AtomicLine2, ...)`
but: the decay rate is `A_rad` times `pesc` (no floor), the oscillator strength is the stored `gf`, the line
energy is that of the wavelength (not of the levels), the photoexcitation is not reduced by a covering
fraction and no line is excluded for its wavelength. `frate` is the decay and `irate` the photoexcitation
from `radiation` at the line energy (zero without it); `init` is the higher and `final` the lower level.
`ienergy` is the photoexcitation rate times the line energy and `fenergy` is 0 (`ucalc` returns only
that one); `opacity` is the line-centre opacity. The level energies only order the two levels. The line
opacity added to the continuum arrays (`linopac`) is not included. The level table needs the Fe UTA levels, which
`level_table` includes.
"""
function rate(coef::RadiativeFeDecay, cell::Cell; mass=atomic_mass(levels_of(coef), coef.ion), vturb=1.0, pesc=1.0, radiation=nothing, index=false, verbose=false)
    levels = levels_of(coef)
    nlev = nlevels(levels, coef.ion)

    K = constants()
    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0., opacity=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 >= nlev || i2 <= 0 || i2 >= nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    up, lo = a.E < b.E ? (b, a) : (a, b)
    index && return (; none..., init=up.level, final=lo.level)

    elin = abs(Float64(coef.λ))
    vtherm = sqrt((vturb*cm_per_km)^2 + (K.thermal_speed/sqrt(mass/cell.T))^2)
    sigma = K.line_xsec*Float64(coef.gf)*elin*cm_per_Å/vtherm
    ener = K.hc_eVÅ/elin
    frate = Float64(coef.A_rad)*pesc
    irate = radiation === nothing ? 0.0 :
        sigma*radiation.F[nbin(radiation, ener)]*vtherm/K.light_speed
    (; init=up.level, final=lo.level, frate, irate, fenergy=0., ienergy=irate*ener*K.ergsev, opacity=sigma)
end

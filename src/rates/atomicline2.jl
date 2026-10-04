# Line (k−i) radiation rates of N-electron ion: r1 = λ(Å); r2 = g f (i,k);
# r3= A(k,i) (s−1); i1= i (lower level); i2= k (upper level); i3= Z; i4= ionN
#
# The stored order of the two levels is not reliable (13% of the records list
# the higher-energy level first), so `rate` orders them by level energy as
# XSTAR's ucalc does.

# XSPEC data type: 50

const AtomicLine2Desc = "op line rad. rates"

const min_wavelength = 1e-34           # Å; lines below this have no wavelength
const no_wavelength = 0.99e9           # Å; longer wavelengths mark lines without opacity
const decay_floor = 1e-20              # floor of the decay rate, per unit ntot
const cm_per_Å = 1e-8
const cm_per_km = 1e5

struct AtomicLine2{I, R, L<:LevelTable} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # the two levels, in stored order
    Z::I      # atomic number
    ion::I    # ion index (XSTAR ionN)
    λ::R      # wavelength (Å)
    gf::R     # weighted oscillator strength
    A::R      # Einstein A (s⁻¹)
    levels::L              # the level data of the database (a LevelTable)
end

function AtomicLine2(rate::Int32, label::String, ivec::I, rvec::R, levels::LevelTable) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AtomicLine2(Int8(rate), label::String, Transition(ivec[1], ivec[2]), ivec[3:4]..., rvec..., levels)
end

"""
    rate(coef::AtomicLine2, cell; mass=atomic_mass(coef.levels, coef.ion), vturb=1.0, pesc=1.0, radiation=nothing, cfrac=0.0, index=false)

Radiative decay and photoexcitation of an atomic line (XSTAR ucalc type 50).
The level data and the number of levels of the ion (the continuum level is excluded) come from the
coefficient (its `levels` field), and so by default does the atomic `mass` (amu) of the element. `vturb` is
the turbulent velocity (km/s), `pesc` the sum of the two escape probabilities (1
when optically thin), `radiation` the incident spectrum (none by default) and `cfrac` the covering
fraction.

Returns `init` (upper) and `final` (lower) levels in energy order, the decay rate
`frate = A·pesc` (s⁻¹, floored at `1e-20 ntot`), the photoexcitation rate `irate`
from `radiation` at the line energy (zero without `radiation`, for lines without a
wavelength, and reduced by `1 - cfrac`), the energies `fenergy` and `ienergy` they
carry (erg s⁻¹ per ion, from the energy difference of the two levels) and the
line-centre `opacity` (zero for lines without a wavelength). Records with no
wavelength (λ = 0) and records with a level missing from the level table give no rates.
The line opacity that XSTAR adds to its continuum arrays (`linopac`) is not
included.
"""
function rate(coef::AtomicLine2, cell::Cell; mass=atomic_mass(coef.levels, coef.ion), vturb=1.0, pesc=1.0, radiation=nothing, cfrac=0.0, index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.,
        opacity=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 >= nlev || i2 <= 0 || i2 >= nlev) && return none
    elin = abs(Float64(coef.λ))
    elin <= min_wavelength && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    up, lo = a.E < b.E ? (b, a) : (a, b)
    index && return (; none..., init=up.level, final=lo.level)

    A = Float64(coef.A)
    frate = max(A*pesc, decay_floor*cell.ntot)
    flin = K.A_to_f*A*up.g*elin^2/lo.g
    vtherm = sqrt((vturb*cm_per_km)^2 + (K.thermal_speed/sqrt(mass/cell.T))^2)
    sigma = K.line_xsec*flin*elin*cm_per_Å/vtherm
    ener = abs(Float64(up.E) - Float64(lo.E))
    irate = 0.0
    if radiation !== nothing && elin <= no_wavelength
        irate = sigma*radiation.F[nbin(radiation, ener)]*vtherm/K.light_speed*max(0.0, 1 - cfrac)
    end
    opacity = elin > no_wavelength ? 0. : sigma
    (; init=up.level, final=lo.level, frate=frate, irate=irate,
        fenergy=frate*ener*K.ergsev, ienergy=irate*ener*K.ergsev, opacity=opacity)
end

# Line (k−i) radiation rates of N-electron ion: r1 = λ(Å); r2 = g f (i,k);
# r3= A(k,i) (s−1); i1= i (lower level); i2= k (upper level); i3= Z; i4= ionN
#
# The stored order of the two levels is not reliable (13% of the records list
# the higher-energy level first), so `rate` orders them by level energy as
# XSTAR's ucalc does.

# XSPEC data type: 50

const AtomicLine2Desc = "op line rad. rates"

struct AtomicLine2{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # the two levels, in stored order
    Z::I      # atomic number
    ion::I    # ion index (XSTAR ionN)
    λ::R      # wavelength (Å)
    gf::R     # weighted oscillator strength
    A::R      # Einstein A (s⁻¹)
end

function AtomicLine2(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AtomicLine2(Int8(rate), label::String, Transition(ivec[1], ivec[2]), ivec[3:4]..., rvec...)
end

"""
    rate(coef::AtomicLine2, cell; levels, mass, vturb=1.0, pesc=1.0, nlev=typemax(Int), index=false)

Radiative decay of an atomic line (XSTAR ucalc type 50). `levels` is a
`level_table`, `mass` the atomic mass (amu) of the element, `vturb` the
turbulent velocity (km/s), `pesc` the sum of the two escape probabilities
(1 when optically thin) and `nlev` the number of levels of the ion (the
continuum level `nlev` is excluded).

Returns `init` (upper) and `final` (lower) levels in energy order, the decay rate
`frate = A·pesc` (s⁻¹, floored at `1e-20 ntot`), the emitted power `fenergy` (erg s⁻¹
per ion, from the energy difference of the two levels) and the line-centre
`opacity` (zero for lines without a wavelength). `irate` and `ienergy` are 0:
XSTAR's photoexcitation from the radiation field is not included until a
radiation object exists. Records with no wavelength (λ = 0) and records with a
level missing from `levels` give no rates.
"""
function rate(coef::AtomicLine2, cell::Cell; levels, mass, vturb=1.0, pesc=1.0,
    nlev=typemax(Int), index=false, verbose=false)

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.,
        opacity=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 >= nlev || i2 <= 0 || i2 >= nlev) && return none
    elin = abs(Float64(coef.λ))
    elin <= 1e-34 && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    up, lo = a.E < b.E ? (b, a) : (a, b)
    index && return (; none..., init=up.level, final=lo.level)

    A = Float64(coef.A)
    frate = max(A*pesc, 1e-20*cell.ntot)
    flin = 1e-16*A*up.g*elin^2/(0.667274*lo.g)
    vtherm = sqrt((vturb*1e5)^2 + (1.29e6/sqrt(mass/cell.T))^2)
    opacity = elin > 0.99e9 ? 0. : 0.02655*flin*elin*1e-8/vtherm
    fenergy = frate*abs(Float64(up.E) - Float64(lo.E))*ergsev
    (; init=up.level, final=lo.level, frate=frate, irate=0., fenergy=fenergy,
        ienergy=0., opacity=opacity)
end

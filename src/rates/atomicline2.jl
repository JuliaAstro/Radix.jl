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
    rate(coef::AtomicLine2, cell; levels, mass, vturb=1.0, index=false)

Radiative decay of an atomic line (XSTAR ucalc type 50). `levels` is a
`level_table`, `mass` the atomic mass (amu) of the element and `vturb` the
turbulent velocity (km/s). Returns `init` (upper) and `final` (lower) levels,
the decay rate `frate` (= A, s⁻¹), the emitted power `fenergy` (erg s⁻¹ per ion)
and the line-centre `opacity` (zero for two-photon records and for lines
without a wavelength). `irate` is always 0: stimulated terms are handled by the
radiation field. If either level is missing from `levels` the rates are 0. For
records without a wavelength (λ = 0) the emitted power uses the level energy
difference, where XSTAR would divide by zero.
"""
function rate(coef::AtomicLine2, cell::Cell; levels, mass, vturb=1.0,
    index=false, verbose=false)

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.,
        opacity=0.)
    a = get(levels, (coef.ion, coef.transition.lower), nothing)
    b = get(levels, (coef.ion, coef.transition.upper), nothing)
    (a === nothing || b === nothing) && return none
    up, lo = a.E < b.E ? (b, a) : (a, b)
    index && return (; none..., init=up.level, final=lo.level)

    A = Float64(coef.A)
    λ = Float64(coef.λ)
    elin = abs(λ)
    flin = 1e-16*A*up.g*elin^2/(0.667274*lo.g)
    vtherm = max(vturb*1e5, 1.3e6/sqrt(mass/cell.T))
    sigma = 0.02655*flin*elin*1e-8/vtherm
    opacity = (λ > 0.99e9 || coef.rtype == 9) ? 0. : sigma
    # photon energy (eV); 339 records have no wavelength (λ = 0), where ucalc
    # would divide by zero, so use the level energies instead
    E = elin > 0 ? 12398.54/elin : up.E - lo.E
    fenergy = A*E*ergsev
    (; init=up.level, final=lo.level, frate=A, irate=0., fenergy=fenergy,
        ienergy=0., opacity=opacity)
end

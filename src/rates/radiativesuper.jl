# Radiative transition rates from superlevels to spectroscopic levels: r1−rjnd
# = (ne( j), j= 1, jnd); rjnd+1−rjnd+nt = (Te( j), j= 1, jnt);
# rjnd+nt+1−rjnd+nt+nt*nd = ((A( j, j′), j′ = 1, j′nd), j= 1, jnt);
# rnd+nt+nt*nd+1 = λ(Å); i1= nd; i2= nt; i3= i (lower level);
# i4= k (upper level); i5= Z; i6= ionN

# XSTAR data type: 71

const RadiativeSuperDesc = "Transition rates from superlevel to spect. lvls"

struct RadiativeSuper{I, R, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    λ::R                # wavelength (Å)
    ne_grid::Vector{R}  # log₁₀ densities (cm⁻³)
    T_grid::Vector{R}   # log₁₀ temperatures (K)
    A::Matrix{R}        # log₁₀ Einstein A (s⁻¹; a single entry above 30 is A itself), size (T, ne)
    levels::L              # the level data of the database (a Levels)
end

function RadiativeSuper(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Nd, Nt = ivec[1:2]
    RadiativeSuper(Int8(rate), label, Transition(ivec[3], ivec[4]), ivec[5:6]..., rvec[Nd+Nt+Nd*Nt+1],
        Vector(rvec[1:Nd]), Vector(rvec[Nd+1:Nd+Nt]),
        reshape(rvec[Nd+Nt+1:Nd+Nt+Nd*Nt], Int(Nt), Int(Nd)), levels)   # ucalc reads the temperature fastest
end

const superdecay_cap = 1e10               # decay rate cap for ions 96 and 97 (calcium I and II)
const superdecay_capped_ions = (96, 97)
const superdecay_A_linear = 30.0          # a one-point table above this holds A itself

# calt71: the Einstein A of the record at the conditions of the cell and its wavelength
function superlevel_decay(coef::RadiativeSuper, T, den)
    logne, logT = Float64.(coef.ne_grid), Float64.(coef.T_grid)
    logA = Float64.(coef.A)
    nden, ntem = length(logne), length(logT)
    if nden == 1 && ntem == 1
        a = logA[1]
        return (10^(a > superdecay_A_linear ? log10(a) : a), Float64(coef.λ))
    end
    (exp10(interpolate_super_table(logne, logT, logA, T, den)), Float64(coef.λ))
end

"""
    rate(coef::RadiativeSuper, cell; vturb=1.0, ptmp=(0.5, 0.5), index=false)

Radiative decay between a superlevel and a spectroscopic level (XSTAR ucalc type 71). The Einstein
A is interpolated in log₁₀ T and log₁₀ n from the table of the record (the temperature may exceed the
table by one dex) and multiplied by the sum `ptmp` of the escape probabilities, with a cap of 10¹⁰ s⁻¹
for ions 96 and 97. As in `ucalc`, the weights come from the first level of the record (`transition.lower`,
`init`) and the decay is returned as `irate` (`frate` is 0). `ienergy` is the decay times the energy of the
line (from its wavelength above 0.1 Å, from the levels otherwise) and `opacity` the line-centre
opacity for the atomic `mass` (amu, from the level table by default) and the turbulent velocity `vturb` (km/s), both zero for lines
longer than 10⁹ Å. Needs levels `init` and `final` between 1 and `nlev`.
"""
function rate(coef::RadiativeSuper, cell::Cell; mass=atomic_mass(coef.levels, coef.ion), vturb=1.0, ptmp=(0.5, 0.5), index=false, kw...)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0., opacity=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    index && return (; none..., init=i1, final=i2)
    (i1 <= 0 || i1 > nlev || i2 <= 0 || i2 > nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    A, elin = superlevel_decay(coef, cell.T*T_unit, cell.ntot)
    decay = A*(ptmp[1] + ptmp[2])
    coef.ion in superdecay_capped_ions && (decay = min(decay, superdecay_cap))
    flin = K.A_to_f*A*a.g*elin^2/b.g
    vtherm = sqrt((vturb*cm_per_km)^2 + (K.thermal_speed/sqrt(mass/cell.T))^2)
    sigma = elin > no_wavelength ? 0. : K.line_xsec*flin*elin*cm_per_Å/vtherm
    ener = abs(Float64(a.E) - Float64(b.E))
    energy = decay*ener*K.ergsev
    elin > 0.1 && (energy = decay*K.hc_eVÅ/(elin + tiny)*K.ergsev_decay)
    (; init=i1, final=i2, frate=0., irate=decay, fenergy=0., ienergy=energy, opacity=sigma)
end

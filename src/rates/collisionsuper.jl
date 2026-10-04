# Collision transition rates from superlevels to spectroscopic levels:
# r1−rjnd = (ne( j), j= 1, jnd); rjnd+1−rjnd+nt = (Te( j), j= 1, jnt);
# rjnd+nt+1−rjnd+nt+nt*nd = ((C( j, j′), j′ = 1, j′ nd), j= 1, jnt) (s−1);
# rjnd+nt+nt*nd+1 = λ(Å); i1= nd; i2= nt; i3= i (lower level);
# i4= k (upper level); i5= Z; i6= ionN

# XSTAR data type: 77

const CollisionSuperDesc = "coll rates from 71"

@with_levels struct CollisionSuper{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    λ::R                # wavelength (Å)
    ne_grid::Vector{R}  # log₁₀ densities (cm⁻³)
    T_grid::Vector{R}   # log₁₀ temperatures (K)
    C::Matrix{R}        # log₁₀ collision rates (s⁻¹), size (T, ne)
end

function CollisionSuper(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Nd, Nt = ivec[1:2]
    CollisionSuper(Int8(rate), label, Transition(ivec[3], ivec[4]), ivec[5:6]..., rvec[Nd+Nt+Nd*Nt+1],
        Vector(rvec[1:Nd]), Vector(rvec[Nd+1:Nd+Nt]),
        reshape(rvec[Nd+Nt+1:Nd+Nt+Nd*Nt], Int(Nt), Int(Nd)))   # ucalc reads the temperature fastest
end

const supercoll_min_dE = 1.0              # eV; no rates for levels closer than this
const supercoll_xt_max = 100.0            # no excitation rate above this ΔE/kT

# the statistical weight 2(2l+1) of a level numbered by shells: k(k−1)/2+1 .. k(k+1)/2 hold l = 0 .. k−1
function supercoll_weight(level)
    k = 0
    while true
        k += 1
        level >= (k + 1)*k ÷ 2 + 1 || break
    end
    2*(2*(level - (k*(k - 1) ÷ 2 + 1)) + 1)
end

"""
    rate(coef::CollisionSuper, cell; index=false)

Collisional transitions between a superlevel and a spectroscopic level (XSTAR ucalc type 77).
`irate` is the de-excitation rate `10^C` interpolated in log₁₀ T and log₁₀ n from the table (the
temperature floored at `2.88e6/λ` K from the level energies), `frate` the excitation rate from its
Boltzmann factor at the wavelength `λ` of the record and the weight `2(2l + 1)` where `l` follows
from numbering the first level of the record by shells. `fenergy` and `ienergy` are the rates times
the energy difference of the levels. Both levels must be in `1:nlev`, differ, and be at least 1 eV apart.
"""
function rate(coef::CollisionSuper, cell::Cell; index=false, verbose=false)
    levels = levels_of(coef)
    nlev = nlevels(levels, coef.ion)
    K = constants()
    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    index && return (; none..., init=i1, final=i2)
    (i1 <= 0 || i1 > nlev || i2 <= 0 || i2 > nlev || i1 == i2) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    dE = Float64(b.E) - Float64(a.E)
    abs(dE) < supercoll_min_dE && return none
    T = max(cell.T*T_unit, K.T_floor_coeff*(dE + tiny)/K.hc_eVÅ_single)
    logC = interpolate_super_table(Float64.(coef.ne_grid), Float64.(coef.T_grid), Float64.(coef.C), T,
        cell.ntot; guard=1e-36)
    cul = exp10(logC)
    xt = K.hc_over_k/Float64(coef.λ)/T
    clu = xt < supercoll_xt_max ? cul*exp(-xt)/supercoll_weight(i1) : zero(cul)
    (; init=i1, final=i2, frate=clu, irate=cul, fenergy=clu*abs(dE)*K.ergsev, ienergy=cul*abs(dE)*K.ergsev)
end

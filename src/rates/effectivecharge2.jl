# Effective charge to be used in collisional ionization, with the hydrogenic formula of Sampson and Zhang (szirco): the record is that of type 57
# r1 = Zeff; i1= n; i2= L; i3= 2J; i4= Z; i5= i; i6= ionN

# XSTAR data type: 65

const EffectiveCharge2Desc = "effective charge to be used in coll. ion."

struct EffectiveCharge2{I, R, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I      # principal quantum number
    L::I      # orbital angular momentum
    twoJ::I   # 2J
    Z::I      # atomic number
    level::I  # level index
    ion::I    # ion index (XSTAR ionN)
    Zeff::R   # effective charge
    levels::L              # the level data of the database (a Levels)
end

function EffectiveCharge2(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    iv = length(ivec) < 6 ? (ivec[1:3]..., Int32(0), ivec[4:5]...) : ivec
    EffectiveCharge2(Int8(rate), label, iv..., rvec[1], levels)
end

const effcharge2_rinf = Float64(2.08f-22)       # the factor of the recombination rate

"""
    rate(coef::EffectiveCharge2, cell; index=false)

Collisional ionization of a level and its inverse, three-body recombination (XSTAR ucalc type 65), from `szirco`, the hydrogenic formula of Sampson and Zhang (1988) for the principal quantum
number `n` and the effective charge `Zeff` of the record. `frate` is the ionization rate (the coefficient times the electron density), `irate` the recombination rate, which follows from it by the Saha factor
(the weights of the ground level and of the continuum, `2.08×10⁻²² g₁/g_c T⁻³ᐟ² e^{ΔE/kT}`), `init` the level of the record and `final` the continuum level (`nlev`). `fenergy` and `ienergy` are the rates
times the ionization potential `ΔE` of the level (ucalc's `ans6` and `ans5`, positive, unlike those of type 57). There is no check of the level: its potential is limited to 0 from below.
"""
function rate(coef::EffectiveCharge2, cell::Cell; index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()
    none = (; init=Int32(0), final=Int32(0), frate=0., irate=0., fenergy=0., ienergy=0.)
    idest1, idest2 = Int(coef.level), nlev
    index && return (; none..., init=idest1, final=idest2)
    lo = get(levels, (coef.ion, idest1), nothing)
    cont = get(levels, (coef.ion, nlev), nothing)
    ground = get(levels, (coef.ion, 1), nothing)
    (lo === nothing || cont === nothing || ground === nothing) && return none
    eth = max(0.0, Float64(cont.E) - Float64(lo.E))
    cii = szirco(Int(coef.n), cell.T*T_unit, Float64(coef.Zeff))
    frate = cii*cell.nₑ
    rinf = effcharge2_rinf*Float64(ground.g)/Float64(cont.g)/cell.T/sqrt(cell.T)
    irate = frate*rinf*expo(eth/(K.kT_eV*cell.T))
    (; init=idest1, final=idest2, frate, irate, fenergy=frate*eth*K.ergsev, ienergy=irate*eth*K.ergsev)
end

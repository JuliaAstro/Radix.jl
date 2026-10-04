# Effective ion charge for i-th level of N-electron ion: r1 = Zeff; i1= n;
# i2= L; i3= 2J; i4= Z; i5= i; i6= ionN

# XSTAR data type: 57

const EffectiveChargeDesc = "effective charge to be used in coll. ion."

struct EffectiveCharge{I, R, L<:LevelTable} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I      # principal quantum number
    L::I      # orbital angular momentum
    twoJ::I   # 2J
    Z::I      # atomic number
    level::I  # level index
    ion::I    # ion index (XSTAR ionN)
    Zeff::R   # effective charge
    levels::L              # the level data of the database (a LevelTable)
end

function EffectiveCharge(rate::Int32, label::String, ivec::I, rvec::R, levels::LevelTable) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    # records of five integers leave out the atomic number
    iv = length(ivec) < 6 ? (ivec[1:3]..., Int32(0), ivec[4:5]...) : ivec
    EffectiveCharge(Int8(rate), label, iv..., rvec[1], levels)
end

const effcharge_n_max_density = 1e18             # cm⁻³
const effcharge_tmin_coeff = 3.8e4
const effcharge_rno1 = 1.8887e8
const effcharge_rno1_exp = Float64(0.3333f0)
const effcharge_rno2 = 1.814e26
const effcharge_rno2_exp = Float64(0.13333f0)
const effcharge_floor = 1e-30
const effcharge_ete_floor = 1e-20
const effcharge_cion_floor = 1e-24
const effcharge_g_floor = 1e-48

# calt57: collisional ionization and three-body recombination rate coefficients of a level of principal quantum number n
function effective_charge_rates(T, den, E, E_ion, n)
    K = constants()
    E_ion < E && return (0.0, 0.0)
    rio = (E_ion - E)/K.Ry_eV_coarse
    rc = sqrt(rio)*n
    den = min(den, effcharge_n_max_density)
    tmin = effcharge_tmin_coeff*rc*sqrt(rc)
    temp = max(T, tmin)
    rno = min(sqrt(effcharge_rno1*rc/den^effcharge_rno1_exp),
        (effcharge_rno2*rc^6/2/den)^effcharge_rno2_exp)
    trunc(Int, rno) > n || return (0.0, 0.0)
    cion = irc(n, temp, rc, rno)
    if T < tmin
        beta = (sqrt((100rc + 91)/(4rc + 3)) - 5)/4
        wte = log(1 + T/K.Ry_K_cb/rio)^(beta/(1 + T/K.Ry_K_cb*rio))
        wtm = log(1 + tmin/K.Ry_K_cb/rio)^(beta/(1 + tmin/K.Ry_K_cb*rio))
        ete = eint(rio/T*K.Ry_K_cb)[1]
        ete < effcharge_ete_floor && return (0.0, 0.0)
        etm = eint(rio/tmin*K.Ry_K_cb)[1]
        cion = cion*sqrt(tmin/T)*ete/(etm + effcharge_floor)*wte/(wtm + effcharge_floor)
    end
    cion <= effcharge_cion_floor && return (cion, 0.0)
    cion /= n*n
    (cion, cion*K.saha_coeff_cb*exp((E_ion - E)*K.eV_K/T)/T^1.5)
end

"""
    rate(coef::EffectiveCharge, cell; index=false)

Collisional ionization of a level and its inverse, three-body recombination (XSTAR ucalc type
57), from the hydrogenic fits of `irc`/`szirc` with the effective charge set by the ionization
potential of the level. `frate` is the ionization rate, `irate` the recombination rate (both
including the electron density), `init` the level and `final` the continuum (`nlev`). `fenergy` and
`ienergy` are the rates times the ionization potential, negated (ucalc's `ans6` and `ans5`). The
ground level, and a level with a non-positive potential, give nothing. The stored `Zeff` is not used.
"""
function rate(coef::EffectiveCharge, cell::Cell; index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()
    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    idest1, idest2 = Int(coef.level), nlev
    index && return (; none..., init=idest1, final=idest2)
    (coef.n <= 0 || idest1 <= 1 || idest1 > nlev) && return none
    lo = get(levels, (coef.ion, idest1), nothing)
    cont = get(levels, (coef.ion, nlev), nothing)
    (lo === nothing || cont === nothing) && return none
    eth = max(0.0, Float64(cont.E) - Float64(lo.E))
    eth <= 0 && return none
    cion, crec = effective_charge_rates(cell.T*T_unit, cell.nₑ, Float64(lo.E), eth, Int(coef.n))
    frate = cion*cell.nₑ
    irate = crec*Float64(lo.g)/(effcharge_g_floor + Float64(cont.g))*cell.nₑ*cell.nₑ
    (; init=idest1, final=idest2, frate, irate, fenergy=-frate*eth*K.ergsev, ienergy=-irate*eth*K.ergsev)
end

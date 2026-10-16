# Electron-impact collision strengths of the SPEX database: the continuum (a formula of `itype` 1 to 8 with six parameters) and the resonances (Gaussian-like terms).
# r1 = the transition energy (keV); r2−r7 = the parameters a₁−a₆ of the continuum; r8−r11 = not used; r12 onwards = the positions x of the nreso resonances (1 to 4, as many as
# `nreso`), then their amplitudes; i1= itype; i2= nreso; i3= i (a level); i4= j (the other level); i5= Z; i6= ionN

# XSTAR data type: 102

const CollisionSPEXDesc = "SPEX collisional data"

struct CollisionSPEX{I, R, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    itype::I            # the formula of the continuum (1 to 8)
    nreso::I            # the number of resonances
    transition::Transition{I}   # the two levels, in stored order
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    E::R                # transition energy (keV)
    a::NTuple{6, R}     # parameters of the continuum
    reso::Vector{R}     # r12 onwards: the positions and amplitudes of the resonances
    levels::L              # the level data of the database (a Levels)
end

function CollisionSPEX(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionSPEX(Int8(rate), label, ivec[1], ivec[2], Transition(ivec[3], ivec[4]), ivec[end - 1], ivec[end], rvec[1],
        (rvec[2], rvec[3], rvec[4], rvec[5], rvec[6], rvec[7]), Vector(rvec[12:end]), levels)
end

# exintn: the exponential integrals Eₙ(x), n = 1 to 6, of the rational approximations of XSTAR (those of Abramowitz and Stegun 5.1.53 and 5.1.56 in single precision) and the recurrence
const exintn_a = Float64.(Float32[7.122452e-07, 1.766345e-06, 2.928433e-05, .0002335379, .001664156, .01041576, .05555682, .2500001, .9999999, .57721566490153])
const exintn_b = Float64.(Float32[8.5733287401, 18.059016973, 8.6347608925, .2677737343, 9.5733223454, 25.6329561486, 21.0996530827, 3.9584969228])

function exintn(x, n)
    x == 0 && return zero(x)
    a, b = exintn_a, exintn_b
    if x <= 1
        e1 = ((((((((a[1]*x - a[2])*x + a[3])*x - a[4])*x + a[5])*x - a[6])*x + a[7])*x - a[8])*x + a[9])*x - log(x) - a[10]
    else
        x2 = x*x
        x3 = x2*x
        x4 = x2*x2
        e1 = (x4 + b[1]*x3 + b[2]*x2 + b[3]*x + b[4])/(x4 + b[5]*x3 + b[6]*x2 + b[7]*x + b[8])
        e1 = e1/(x*expo(x))
    end
    n == 1 && return e1
    n == 2 && return expo(-x) - x*e1
    n == 3 && return (x*x*e1 + expo(-x)*(1 - x))/2
    n == 4 && return (expo(-x)*(2 - x + x*x) - x*x*x*e1)/6
    n == 5 && return (expo(-x)*(6 - 2x + x*x - x*x*x) + x*x*x*x*e1)/24
    (expo(-x)*(24 - 6x + 2x*x - x*x*x + x*x*x*x) - x*x*x*x*x*e1)/120
end

# the parameters that each formula (itype 1 to 8) uses, in the order of par(1) to par(9) of `omega_c`; 0: the parameter is 0
const spex_parameters = ((1, 2, 3, 4, 0, 0, 0, 5, 0), (1, 2, 3, 4, 0, 0, 0, 0, 0), (0, 0, 1, 2, 3, 0, 0, 0, 4), (1, 2, 3, 4, 0, 0, 0, 5, 6),
    (1, 2, 3, 4, 0, 0, 0, 0, 5), (0, 0, 1, 2, 3, 4, 0, 0, 5), (0, 0, 0, 1, 2, 3, 4, 0, 5), (6, 0, 1, 2, 3, 4, 0, 0, 5))

# omega_c: the continuum collision strength (upsilon) at x = ΔE/kT for the formula `itype` with the parameters `a`
function spex_omega_c(x, itype, a)
    ind = spex_parameters[itype]
    par = ntuple(i -> ind[i] == 0 ? 0.0 : Float64(a[ind[i]]), 9)
    xr = x*(1 + par[9])
    ex = exp(xr)
    en = ntuple(n -> ex*exintn(xr, n), 6)
    par[1] + xr*(par[2]*en[1] + par[3]*en[2] + 2*par[4]*en[3] + 6*par[5]*en[4] + 24*par[6]*en[5] + 120*par[7]*en[6]) + par[8]*en[1]
end

# omega_r: the resonance part, a sum of A x e^{-x}
spex_omega_r(x, xreso, areso, nreso) = sum(areso[i]*(xreso[i]*x)*exp(-xreso[i]*x) for i in 1:clamp(nreso, 0, 4); init=0.0)

const spex_min_energy = 1e-34          # keV is multiplied by 1000 and added to this before the wavelength is made
const spex_min_wavelength = 1e-24

"""
    rate(coef::CollisionSPEX, cell; index=false)

Electron-impact excitation and de-excitation from the SPEX collision strengths (XSTAR ucalc type 102): Υ is the sum of the continuum, one of eight formulas (`itype`) of six parameters
in `x = ΔE/kT` with the exponential integrals `E₁` to `E₆`, and of the resonances, `Σ Aᵢ xᵢ x e^{-xᵢ x}` (`nreso`, at most four: ucalc reads the parameters of four whenever `nreso` is not 1, 2 or 3 but sums over `nreso` of them, so that none are added for `nreso` = 0); ΔE is 1000 times the first real (keV). `irate` is the de-excitation rate
`nₑ C Υ/(√T g_u)` and `frate` the excitation rate `irate g_u e^{-ΔE/kT}/g_l` (T in 10⁴ K); `fenergy` and `ienergy` are the rates times ΔE. As in ucalc, `init` is the level of the two with the higher energy, `final` the other one (so that
the rates go the opposite way from those of the other collision types), and nothing limits Υ to positive values. The records need `itype` between 1 and 8.
"""
function rate(coef::CollisionSPEX, cell::Cell; index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()
    none = (; init=Int32(0), final=Int32(0), frate=0., irate=0., fenergy=0., ienergy=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 > nlev || i2 <= 0 || i2 > nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    up, lo = a.E < b.E ? (b, a) : (a, b)
    index && return (; none..., init=up.level, final=lo.level)
    eij = 1000*Float64(coef.E)
    elin = K.hc_eVÅ_single/(spex_min_energy + eij)
    elin <= spex_min_wavelength && return none
    1 <= coef.itype <= 8 || return none
    delt = K.hc_eVÅ_single/(spex_min_energy + elin)/(K.kT_eV*cell.T)
    nreso = Int(coef.nreso)
    n = nreso in (1, 2, 3) ? nreso : 4          # (ucalc reads four resonances for any other number)
    r = coef.reso
    xreso = ntuple(i -> i <= n ? Float64(r[i]) : 0.0, 4)
    areso = ntuple(i -> i <= n ? Float64(r[n + i]) : 0.0, 4)
    omega_cont = spex_omega_c(delt, Int(coef.itype), coef.a)
    omega_res = spex_omega_r(delt, xreso, areso, nreso)
    ggup, gglo = Float64(up.g), Float64(lo.g)
    cji_cont = K.collision_rate_coeff*omega_cont/sqrt(cell.T)/ggup
    cji_res = K.collision_rate_coeff*omega_res/sqrt(cell.T)/ggup
    exptmp = expo(-delt)
    cji = cji_cont + cji_res
    cij = cji_cont*ggup*exptmp/gglo + cji_res*ggup*exptmp/gglo
    frate, irate = cij*cell.nₑ, cji*cell.nₑ
    (; init=up.level, final=lo.level, frate, irate, fenergy=frate*eij*K.ergsev, ienergy=irate*eij*K.ergsev)
end

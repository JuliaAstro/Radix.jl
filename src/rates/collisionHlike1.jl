# Analytic fits for effective collision strengths in H-like ions [14]:
# r1−rjmax = coefficients; i1= i (lower level); i2= k (upper level); i3= 1;
# i8= ionN ; s1= Transition

# XSTAR data type: 60

const CollisionHlike1Desc = ""

struct CollisionHlike1{I, R, L<:LevelTable} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    ion::I             # ion index (XSTAR ionN)
    coeffs::Vector{R}  # fit coefficients
    levels::L              # the level data of the database (a LevelTable)
end

function CollisionHlike1(rate::Int32, label::String, ivec::I, rvec::R, levels::LevelTable) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionHlike1(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[4], Vector(rvec), levels)
end

const hlike_T_high = 1e9          # K; above this the scaled temperature is fixed
const hlike_t_max = 1.0           # the fit is used up to this scaled temperature (kT in Rydberg)
const hlike_g_guard = 1e-16       # added to the weights
const hlike_min_dE = 1e-24

# calt6062: effective collision strength from the polynomial coefficients c = r1, r2, r3, ...: c[3:m] for the
# first form (type 60) and c[3:m-3] plus a log-exponential term for the second (type 62), with the scaled
# temperature kT in Rydberg limited to 1 and the extrapolation of the fit beyond it
function hlike_upsilon(c, T, second_form)
    K = constants()
    m = length(c)
    t1 = T > hlike_T_high ? hlike_T_high*K.Ry_per_K : T*K.Ry_per_K
    tt = min(t1, hlike_t_max)
    if second_form
        ups = sum(c[i + 2]*tt^(i - 1) for i in 1:m - 5; init=zero(tt)) +
            c[m - 2]*log(c[m - 1]*tt)*exp(-c[m]*tt)
    else
        ups = sum(c[i + 2]*tt^(i - 1) for i in 1:m - 2; init=zero(tt))
    end
    t1 > tt ? ups*(1 + log(t1/hlike_t_max)/(log(t1/hlike_t_max) + 1)) : ups
end

# the collision rates of types 60 and 62
function hlike_collision(coef, cell::Cell, index, second_form)
    K = constants()
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 > nlev || i2 <= 0 || i2 > nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    lo, up = b.E < a.E ? (b, a) : (a, b)
    index && return (; none..., init=lo.level, final=up.level)
    dE = abs(Float64(up.E) - Float64(lo.E))
    dE <= hlike_min_dE && return none
    T = max(cell.T*T_unit, dE*T_unit/K.kT_eV/T_floor_dE_over_kT)
    cijpp = hlike_upsilon(Float64.(coef.coeffs), T, second_form)
    gu, gl = Float64(up.g), Float64(lo.g)
    cji = K.collision_rate_coeff*cijpp/sqrt(cell.T)/(hlike_g_guard + gu)
    cij = cji*gu*expo(-dE/(K.kT_eV*cell.T))/(hlike_g_guard + gl)
    frate, irate = cij*cell.nₑ, cji*cell.nₑ
    (; init=lo.level, final=up.level, frate, irate, fenergy=frate*dE*K.ergsev, ienergy=irate*dE*K.ergsev)
end

"""
    rate(coef::CollisionHlike1, cell; index=false)

Collisional excitation and de-excitation of an H-like ion from a polynomial fit of the effective collision
strength in the scaled temperature kT/Ry (XSTAR ucalc type 60, `calt6062`): Υ = Σ cᵢ τ^{i-1} with the
coefficients from the third real on, τ limited to 1 and Υ extrapolated beyond it by `1 + ln τ₁/(ln τ₁ + 1)`.
`irate` is the de-excitation rate `nₑ C Υ/(√T g_u)` and `frate = irate g_u e^{-ΔE/kT}/g_l` (T in 10⁴ K), `fenergy`
and `ienergy` the rates times ΔE. `init` is the lower and `final` the upper level by energy, both in `1:nlev`;
the temperature of the fit is at least ΔE/50k. Weights have 1e-16 added.
"""
rate(coef::CollisionHlike1, cell::Cell; index=false) =
    hlike_collision(coef, cell, index, false)

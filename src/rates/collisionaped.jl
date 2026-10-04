# APED collision strengths [20]: r1 = T_min, r2 = T_max (K); r3−r22 = (Te( j), j= 1, 20) (K);
# r23−r42 = (Υ( j), j= 1, 20); i1= i (lower level); i2= k (upper level); i3 = kind of fit (100 + n:
# n points of an interpolated electron collision strength); i4= Z; i5= ionN

# XSTAR data type: 92

const CollisionAPEDDesc = "aped collision strengths"

struct CollisionAPED{I, R, L<:LevelTable} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    kind::I             # kind of fit (XSTAR coll_type)
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    T_min::R            # lowest temperature of the fit (K)
    T_max::R            # highest temperature of the fit (K)
    T_grid::Vector{R}   # temperatures (K)
    Υ::Vector{R}        # effective collision strengths
    levels::L              # the level data of the database (a LevelTable)
end

function CollisionAPED(rate::Int32, label::String, ivec::I, rvec::R, levels::LevelTable) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    n = (length(rvec) - 2) ÷ 2
    CollisionAPED(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3], ivec[4:5]..., rvec[1], rvec[2],
        Vector(rvec[3:2 + n]), Vector(rvec[3 + n:2 + 2n]), levels)
end

const aped_interpolated_electron = 100            # kinds 100 + n: n points, electron collision strength
const aped_max_points = 20
const aped_max_chi = 200.0                        # no excitation rate above this ΔE/kT
const aped_min_dE = 1e-34
const aped_keV = 1e3

# interpol_huntd: interpolation of y(x) at z, in log-log where both are positive and linearly otherwise,
# and 0 outside the first n − 1 points (no extrapolation; the last one is not used)
function interpolate_aped(n, x, y, z)
    inc = x[n - 1] > x[1]
    (inc ? (z < x[1] || z > x[n - 1]) : (z > x[1] || z < x[n - 1])) && return zero(z*y[1])
    jl, ju = 0, n
    while ju - jl > 1
        jm = (ju + jl) ÷ 2
        (z > x[jm]) == inc ? (jl = jm) : (ju = jm)
    end
    (jl <= 0 || ju <= 0) && return zero(z*y[1])
    if x[jl] > 0 && x[ju] > 0 && y[jl] > 0 && y[ju] > 0
        grad = (log10(y[ju]) - log10(y[jl]))/(log10(x[ju]) - log10(x[jl]))
        10^(log10(y[jl]) + grad*(log10(z) - log10(x[jl])))
    else
        y[jl] + (y[ju] - y[jl])*((z - x[jl])/(x[ju] - x[jl]))
    end
end

"""
    rate(coef::CollisionAPED, cell; index=false)

Collisional excitation and de-excitation between two levels from a table of the effective collision
strength Υ(T) (XSTAR ucalc type 92, `calc_maxwell_rates`). Only the kinds `100 + n` that the database has
are supported: Υ is interpolated in log-log in `T_grid` (the first `n` points, no extrapolation, zero
outside `T_min`..`T_max` or the grid) and the rates are `frate = 8.629e-6 Υ e^{-χ}/(√T g_l)`,
`irate = 8.629e-6 Υ/(√T g_u)` with χ = ΔE/kT and the electron density; `fenergy` and `ienergy` are the
rates times ΔE. The first level of the record is the lower one whatever the energies. Other kinds throw an error.
"""
function rate(coef::CollisionAPED, cell::Cell; index=false)
    levels = coef.levels
    K = constants()
    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    index && return (; none..., init=i1, final=i2)
    n = Int(coef.kind) - aped_interpolated_electron
    0 < n <= aped_max_points || error("CollisionAPED kind $(coef.kind) is not supported")
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    eij = abs(Float64(b.E) - Float64(a.E))
    eij/aped_keV < aped_min_dE && return none
    T = cell.T*T_unit
    (T < coef.T_min || T > coef.T_max) && return none
    χ = eij/aped_keV/(K.kB_keV*T)
    υ = max(interpolate_aped(n, Float64.(coef.T_grid), Float64.(coef.Υ), T), 0.0)
    frate = irate = zero(υ)
    if χ < aped_max_chi
        frate = K.ups_coeff*υ*exp(-χ)/(sqrt(T)*Float64(a.g))*cell.nₑ
        irate = K.ups_coeff*υ/(sqrt(T)*Float64(b.g))*cell.nₑ
    end
    (; init=i1, final=i2, frate, irate, fenergy=frate*eij*K.ergsev, ienergy=irate*eij*K.ergsev)
end

# Fits to LS collision strengths for He-like ions [32]: r1−rjmax = 
# coefficients; i1= i (lower level); i2= k (upper level); i3= Z; i8= ionN

# XSTAR data type: 69

const CollisionLSDesc = "Kato & Nakazaki (1996) fit to Helike coll. strgt"

struct CollisionLS{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I               # atomic number
    ion::I             # ion index (XSTAR ionN)
    coeffs::Vector{R}  # fit coefficients
end

function CollisionLS(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionLS(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., Vector(rvec))
end

const ls_y_max = 77.0           # the argument of the exponential integral is limited to this
const ls_y_min = 5e-2           # ... and, in type 69, from below to this
const ls_y_huge = 1e20          # no collision strength above this scaled energy
const ls_nine = 9               # number of coefficients of the second form of type 69

# one term of the fits of calt66 and calt69: the scaled energy y = E/kT is clamped to (y_min, 77)
function ls_term(E, a, b, c, d, e, T, y_min)
    K = constants()
    y = clamp(E/T*K.eV_K_ls, y_min, ls_y_max)
    em1 = expint(y)
    y*((a/y + c) + d*0.5*(1 - y)) + em1*(b - c*y + d*y*y*0.5 + e/y)
end

# calt69: Kato and Nakazaki fit; the first form (6 coefficients) or, with 9, a second one with a
# non-resonant part to the threshold x1
function ls_gamma(r, T)
    K = constants()
    y = r[1]/T*K.eV_K_ls
    y > ls_y_huge && return zero(y)
    length(r) == ls_nine || return ls_term(r[1], r[2], r[3], r[4], r[5], r[6], T, ls_y_min)
    y = clamp(y, ls_y_min, ls_y_max)
    a, b, c, d, e, p, q, x1 = r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9]
    em1 = expint(y*x1)
    gnr = a/y + c/x1 + d*0.5*(1/x1^2 - y/x1) + e/y*log(x1) + em1/y/x1*(b - c*y + d*y*y*0.5 + e/y)
    gnr *= y*exp(y*(1 - x1))
    gr = p*(1 + 1/y)*(1 - exp(y*(1 - x1))*(x1 + 1/y)/(1 + 1/y)) + q*(1 - exp(y*(1 - x1)))
    gnr + gr
end

# calt66: fine-structure fit, one term of calt69's first form and, with more than 6 coefficients, two more
function ls_fine_gamma(r, T)
    K = constants()
    y = r[1]/T*K.eV_K_ls
    y > ls_y_huge && return zero(y)
    gam = ls_term(r[1], r[2], r[3], r[4], r[5], r[6], T, 0.0)
    length(r) > 6 || return gam
    gam + ls_term(r[7], r[8], r[9], r[10], r[11], r[12], T, 0.0) +
        ls_term(r[13], r[14], r[15], r[16], r[17], r[18], T, 0.0)
end

# the collision rates of types 66, 68 and 69; `kind` is :fine, :helike or :ls
function helike_fit_collision(coef, cell::Cell, levels, nlev, index, kind)
    K = constants()
    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 > nlev || i2 <= 0 || i2 > nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    lo, up = b.E < a.E ? (b, a) : (a, b)
    index && return (; none..., init=lo.level, final=up.level)
    r = Float64.(coef.coeffs)
    T = cell.T*T_unit
    if kind == :ls
        E = abs(Float64(up.E) - Float64(lo.E))
        Υ = max(ls_gamma(r, T), 0.0)
    elseif kind == :helike
        E = abs(Float64(up.E) - Float64(lo.E) + tiny)
        T = max(T, K.T_floor_coeff*E/K.hc_eVÅ_single)
        tt = log10(T/Float64(coef.Z)^3)
        Υ = max(r[1] + r[2]*tt + r[3]*tt*tt, 0.0)
    else
        E = r[1]                                        # the energy of the record, not of the levels
        E <= tiny && return none
        T = max(T, K.T_floor_coeff*E/K.hc_eVÅ_single)
        Υ = ls_fine_gamma(r, T)
    end
    gu, gl = Float64(up.g), Float64(lo.g)
    cji = K.collision_rate_coeff*Υ/sqrt(cell.T)/gu
    cij = cji*gu*expo(-E/(K.kT_eV*cell.T))/gl
    frate, irate = cij*cell.nₑ, cji*cell.nₑ
    (; init=lo.level, final=up.level, frate, irate, fenergy=frate*E*K.ergsev, ienergy=irate*E*K.ergsev)
end

"""
    rate(coef::CollisionLS, cell; levels, nlev=typemax(Int), index=false)

Collisional excitation and de-excitation of a He-like ion from the fit of Kato and Nakazaki to the
effective collision strength (XSTAR ucalc type 69, `calt69`). The first coefficient is the energy (eV) of the
fit, which scales the temperature, and the next five (or eight, with a non-resonant part) the coefficients; Υ is
clamped at 0 and the temperature is not floored. `irate` is the de-excitation rate `nₑ C Υ/(√T g_u)`,
`frate = irate g_u e^{-ΔE/kT}/g_l` (T in 10⁴ K) and `fenergy`, `ienergy` the rates times the energy difference of the
levels; `init` is the lower and `final` the upper level by energy, both in `1:nlev`.
"""
rate(coef::CollisionLS, cell::Cell; levels, nlev=typemax(Int), index=false, verbose=false) =
    helike_fit_collision(coef, cell, levels, nlev, index, :ls)

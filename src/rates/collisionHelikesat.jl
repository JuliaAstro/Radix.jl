# Fit to effective collision strengths for satellite levels of He-like ions
# [47]: r1−rj7 = fit coefficients; i1= i (lower level); i2= j (upper level);
# i3= Z; i4= ionN

# XSTAR data type: 73

const CollisionHelikeSatDesc =  "Fit to coll. strengths satellite lvls Helike ion"

struct CollisionHelikeSat{I, R, L<:LevelTable} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I               # atomic number
    ion::I             # ion index (XSTAR ionN)
    coeffs::Vector{R}  # fit coefficients
    levels::L              # the level data of the database (a LevelTable)
end

function CollisionHelikeSat(rate::Int32, label::String, ivec::I, rvec::R, levels::LevelTable) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionHelikeSat(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., Vector(rvec), levels)   
end

const sat_gam_high, sat_gam_mid, sat_gam_low = -0.2, 0.0, 0.2   # shifts of the effective charge by 2s weight
const sat_z2s_high, sat_z2s_low = 0.1, 0.01
const sat_y_max = 40.0             # no rate above this Z² E/kT
const sat_exp_max = 80.0           # the exponential integrals are evaluated below this argument (they enter the fit below 40)
const sat_e1_coeff = 1.55

# calt73: Maxwellian collision rate coefficient (cm³ s⁻¹) of the fit at the temperature T (K)
function helike_satellite_rate(Z, c, T)
    K = constants()
    E, z2s, a, co, cr, cr1, r = c
    y = Z^2*E*K.Ry_K_sat/T
    y > sat_y_max && return zero(y)
    gam = z2s >= sat_z2s_high ? sat_gam_high : z2s > sat_z2s_low ? sat_gam_mid : sat_gam_low
    zeff = Z - gam
    e1 = expint(y)/y*exp(-y)
    ya = y*a + y
    ee1, ee2, ee3 = ya <= sat_exp_max ? eint(ya) : (zero(ya), zero(ya), zero(ya))
    er, er1 = r == 1 ? (ee1, ee2) : r == 2 ? (ee2, ee3) : (zero(ya), zero(ya))
    qij = co*exp(-y) + sat_e1_coeff*z2s*e1
    ya <= sat_y_max && (qij += y*exp(y*a)*(cr*er/(a + 1)^(r - 1) + cr1*er1/(a + 1)^r))
    max(qij*K.Ry_K_sat/T*sqrt(T)/zeff^2*K.thermal_bohr, zero(qij))
end

"""
    rate(coef::CollisionHelikeSat, cell; index=false)

Collisional excitation and de-excitation involving a satellite level of a He-like ion from a fit to
the effective collision strength (XSTAR ucalc type 73, `calt73`). The first coefficient is used both as an
energy in Rydberg (in the fit) and as a wavelength in Å (for the temperature floor ΔE/50k, the Boltzmann
factor and the energy of the rates), as `ucalc` does. `irate` is the de-excitation rate, `frate = irate g_u
e^{-ΔE/kT}/g_l`, `fenergy` and `ienergy` the rates times that energy; `init` is the lower and `final` the
upper level by energy, both in `1:nlev`.
"""
function rate(coef::CollisionHelikeSat, cell::Cell; index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()
    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 > nlev || i2 <= 0 || i2 > nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    lo, up = b.E < a.E ? (b, a) : (a, b)
    index && return (; none..., init=lo.level, final=up.level)
    elin = abs(Float64(coef.coeffs[1]))
    elin <= tiny && return none
    dE = K.hc_eVÅ_single/elin                        # (the wavelength reading of the energy)
    T = max(cell.T*T_unit, K.T_floor_coeff/elin)
    crate = helike_satellite_rate(Int(coef.Z), Float64.(coef.coeffs), T)
    cijpp = max(crate/Float64(lo.g), 0.0)
    cji = K.collision_rate_coeff*cijpp/sqrt(cell.T)/Float64(up.g)
    cij = cji*Float64(up.g)*expo(-dE/(K.kT_eV*cell.T))/Float64(lo.g)
    frate, irate = cij*cell.nₑ, cji*cell.nₑ
    (; init=lo.level, final=up.level, frate, irate, fenergy=frate*dE*K.ergsev, ienergy=irate*dE*K.ergsev)
end

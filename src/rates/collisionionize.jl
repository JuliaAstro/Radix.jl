# Collisional ionization rates for N-electron ion [11]: r1 = E(th) (eV);
# r2 = T0 (K); then the scaled temperature grid x(j) = 1 − ln 2/ln(kT/E(th) + 2), j = 1..jmax, and the
# effective collision strengths ρ(j), j = 1..jmax (the header of the source omits the grid);
# i1= i (level); i2 = ionN

# XSTAR data type: 95

const CollisionIonizeDesc = "Bryans CI rates"

struct CollisionIonize{I, R, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    level::I      # level index
    ion::I        # ion index (XSTAR ionN)
    E_th::R       # threshold energy (eV)
    T0::R         # K
    x_grid::Vector{R}  # scaled temperature of the fit
    ρ::Vector{R}  # effective collision strengths
    levels::L              # the level data of the database (a Levels)
end

function CollisionIonize(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    n = (length(rvec) - 2) ÷ 2
    CollisionIonize(Int8(rate), label, ivec..., rvec[1], rvec[2],
        Vector(rvec[3:2 + n]), Vector(rvec[3 + n:2 + 2n]), levels)
end

const ci_ln2 = Float64(0.693147f0)    # ln 2 as ucalc writes it (single precision)
const ci_units = 1e-6                 # the rate is 1e-6 E₁ ρ / √(τ E³) (cm³ s⁻¹ with E in eV)

"""
    rate(coef::CollisionIonize, cell; index=false)

Collisional ionization of an ion and its inverse, three-body recombination (XSTAR ucalc type 95), from the
effective collision strength ρ(x) of the record, interpolated linearly in the scaled temperature
`x = 1 - ln 2/ln(kT/E_th + 2)`. `frate` and `irate` include the electron density; `init` is the ground
level (`ucalc` takes the statistical weight and the level of the ground state whatever the level of the
record) and `final` the continuum, `nlev` (records of rate type 15 give level 1 again). `fenergy` and `ienergy` are the
rates times the threshold energy. For temperatures below the first point of the table ucalc interpolates
with the first interval's left end read from `T0` and the last `x_grid`; reproduced.
"""
function rate(coef::CollisionIonize, cell::Cell; index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()
    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    final = coef.rtype == 5 ? nlev : 1
    index && return (; none..., init=1, final)
    ground = get(levels, (coef.ion, 1), nothing)
    cont = get(levels, (coef.ion, nlev), nothing)
    (ground === nothing || cont === nothing) && return none
    ee = Float64(coef.E_th)
    tt = K.kT_eV*cell.T/ee
    xx = 1 - ci_ln2/log(tt + 2)
    n = length(coef.x_grid)
    # x(0) and y(0) are read from the record before the table: T0 and the last x
    x(m) = m == 0 ? Float64(coef.T0) : Float64(coef.x_grid[m])
    y(m) = m == 0 ? Float64(coef.x_grid[n]) : Float64(coef.ρ[m])
    m = 1
    while m < n && xx > x(m)
        m += 1
    end
    ρ = y(m - 1) + (xx - x(m - 1))*(y(m) - y(m - 1))/(x(m) - x(m - 1))
    e1 = eint(1/tt)[1]
    frate = ci_units*e1*ρ/sqrt(tt*ee^3)*cell.nₑ
    irate = frate*K.saha_ci*Float64(ground.g)/Float64(cont.g)/cell.T/sqrt(cell.T)*cell.nₑ/expo(-1/tt)
    (; init=1, final, frate, irate, fenergy=frate*ee*K.ergsev, ienergy=irate*ee*K.ergsev)
end

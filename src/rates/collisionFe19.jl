# Collision strengths for Fe XIX [10]: r1 = Υ(k,i); i1= i (lower level);
# i2= k (upper level); i3= Z; i4= ionN

# XSTAR data type: 81

const CollisionFe19Desc = "Bhatia Fe XIX collision strengths"

struct CollisionFe19{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I        # atomic number
    ion::I      # ion index (XSTAR ionN)
    Υ::R  # effective collision strength
end

function CollisionFe19(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionFe19(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., rvec[1])    
end

"""
    rate(coef::CollisionFe19, cell; levels, nlev=typemax(Int), index=false)

Collisional excitation and de-excitation of Fe XIX from a constant effective collision strength
Υ (XSTAR ucalc type 81): `irate = nₑ C Υ/(√T g_u)` and `frate = irate g_u e^{-ΔE/kT}/g_l` (T in 10⁴ K, `C` the
Maxwellian coefficient), `fenergy` and `ienergy` the rates times ΔE. `init` is the lower and `final`
the upper level by energy (the stored order is not reliable); both must be in `1:nlev`. Negative Υ
is taken as 0.
"""
function rate(coef::CollisionFe19, cell::Cell; levels, nlev=typemax(Int), index=false, verbose=false)
    K = constants()
    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 > nlev || i2 <= 0 || i2 > nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    lo, up = b.E < a.E ? (b, a) : (a, b)
    index && return (; none..., init=lo.level, final=up.level)
    ΔE = abs(Float64(up.E) - Float64(lo.E))
    cji = K.collision_rate_coeff*max(Float64(coef.Υ), 0.0)/sqrt(cell.T)/Float64(up.g)
    cij = cji*Float64(up.g)*expo(-ΔE/(K.kT_eV*cell.T))/Float64(lo.g)
    frate, irate = cij*cell.nₑ, cji*cell.nₑ
    (; init=lo.level, final=up.level, frate, irate,
        fenergy=frate*ΔE*K.ergsev, ienergy=irate*ΔE*K.ergsev)
end

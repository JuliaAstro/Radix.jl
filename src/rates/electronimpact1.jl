# Electron-impact effective collision strengths for the k−i transition of
# N-electron ion: r1−rjmax = (log Te( j), j= 1, jmax) (K); rj(max+1)−rj(2*max)=
# (Υ(Te( j)), j= 1, jmax) (effective collision strength); i1= i (lower level);
# i2= k (upper level); i3= Z; i5= ionN

# XSTAR data type: 56

const ElectronImpact1Desc = "tabulated collision strength, bautista"
const min_dE = 1e-16                  # eV; degenerate levels are skipped
const min_upsilon = 1e-48              # floor of the tabulated Υ before interpolating

@with_levels struct ElectronImpact1{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    T_grid::Vector{R}   # log₁₀ of the temperatures (K)
    Υ::Vector{R}       # effective collision strengths
end

function ElectronImpact1(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ElectronImpact1(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., Vector(rvec[1:end÷2]),
        Vector(rvec[end÷2+1:end]))
end

"""
    rate(coef::ElectronImpact1, cell; index=false)

Electron-impact excitation and de-excitation (XSTAR ucalc type 56). The levels come from the
coefficient's level table (`attach_levels`). The levels are ordered by energy, `init` is the lower and
`final` the upper level, `frate` the excitation rate and `irate` the
de-excitation rate (s⁻¹), related by detailed balance. Υ is interpolated
linearly in log T; outside the table the nearest segment is extrapolated (as
ucalc does, since its `hunt` clamps the index) and the result clamped at 0. Both
rates are 0 if a level is missing from the level table or the two energies coincide.
"""
function rate(coef::ElectronImpact1, cell::Cell; index=false, verbose=false)
    levels = levels_of(coef)
    K = constants()

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.)
    a = get(levels, (coef.ion, coef.transition.lower), nothing)
    b = get(levels, (coef.ion, coef.transition.upper), nothing)
    (a === nothing || b === nothing) && return none
    lo, up = b.E < a.E ? (b, a) : (a, b)
    ΔE = abs(Float64(up.E) - Float64(lo.E))
    ΔE <= min_dE && return none
    index && return (; none..., init=lo.level, final=up.level)

    logT = log10(cell.T*T_unit)
    T, Υ = coef.T_grid, coef.Υ
    j = clamp(searchsortedlast(T, logT), 1, max(length(T) - 1, 1))
    Υ0 = max(min_upsilon, Float64(Υ[j]))
    cijpp = length(T) == 1 ? Υ0 :
        (Float64(Υ[j+1]) - Υ0)*(logT - T[j])/(T[j+1] - T[j] + tiny) + Υ0
    cijpp = max(0., cijpp)
    cij = K.collision_rate_coeff*cijpp*expo(-ΔE/(K.kT_eV*cell.T))/sqrt(cell.T)/lo.g
    cji = K.collision_rate_coeff*cijpp/sqrt(cell.T)/up.g
    frate, irate = cij*cell.nₑ, cji*cell.nₑ
    (; init=lo.level, final=up.level, frate, irate,
        fenergy=frate*ΔE*K.ergsev, ienergy=irate*ΔE*K.ergsev)
end

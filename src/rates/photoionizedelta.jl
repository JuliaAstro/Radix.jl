# Delta functions to add to photoionization cross sections to match ADF DR
# rates: r1 = E(∞) (eV); r1−rjm = (E( j), j= 1, jm) (eV); rjm+1−rj2m = 
# ( f ( j), j= 1, jm) (cm2); i1= n; i2= L; i3= 2S + 1; i4= Z; i5= kN−1;
# i6= ionN−1; i7= iN ; i8= ionN

# XSTAR data type: 74

const PhotoionizeDeltaDesc = "Delta functions to add to phot. x-sections  DR"

struct PhotoionizeDelta{I, R, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I               # principal quantum number
    L::I               # orbital angular momentum
    spin_mult::I       # 2S+1
    Z::I               # atomic number
    parent::Parent{I}           # parent ion and level
    level::I           # level index
    ion::I             # ion index (XSTAR ionN)
    E_inf::R           # eV
    E_grid::Vector{R}  # energies (eV)
    f::Vector{R}       # cm²
    levels::L              # the level data of the database (a Levels)
end

function PhotoionizeDelta(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    m = (length(rvec) - 1) ÷ 2
    PhotoionizeDelta(Int8(rate), label, ivec[1:4]..., Parent(ivec[6], ivec[5]),
        ivec[7:8]..., rvec[1], Vector(rvec[2:m + 1]), Vector(rvec[m + 2:2m + 1]), levels)
end

const delta_kB = 1.38066e-16               # erg K⁻¹ (the value used in calt74)
const delta_ryd_K = 4.589343e10            # scale of the exponent of the recombination sum
const delta_exp_limit = 40.0
const delta_alpha_coeff = 213.9577e-9
const delta_rate_coeff = 4.752e-22

"""
    rate(coef::PhotoionizeDelta, cell; nlev)

Delta-function photoionization added to the cross sections to match dielectronic
recombination rates (XSTAR ucalc type 74, `calt74`). `frate` sums the spectrum of
`radiation` at the line energies weighted by the line strengths, `irate` is the matching
recombination sum at the temperature of `cell`, scaled by the statistical weights of the
level and of the continuum level `nlev`. `final` is `nlev`. Both are 0 if the spectrum
does not reach the highest line.
"""
function rate(coef::PhotoionizeDelta, cell::Cell; radiation=NO_RADIATION, index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()

    none = (; init=0, final=0, frate=0., irate=0.)
    idest1 = Int(coef.level)
    index && return (; none..., init=idest1, final=nlev)
    lo = get(levels, (coef.ion, idest1), nothing)
    cont = get(levels, (coef.ion, nlev), nothing)
    (lo === nothing || cont === nothing) && return none
    xt = Float64(coef.E_inf)
    x = Float64.(coef.E_grid)
    hgh = Float64.(coef.f)
    m = length(x)
    m == 0 && return none

    te = cell.T*T_unit*delta_kB
    alpha = 0.0
    for i in 1:m
        x[i]/delta_ryd_K/te < delta_exp_limit &&
            (alpha += exp(-x[i]/delta_ryd_K/te)*(x[i] + xt)^2*hgh[i])
    end
    alpha *= delta_alpha_coeff/te^1.5/delta_ryd_K^2

    E, F = radiation.E, radiation.F
    np = length(E)
    E[np] < (x[m] + xt)*K.Ry_eV && return (; init=idest1, final=nlev, frate=0., irate=0.)
    # the spectrum at the line energies, interpolated linearly in the grid
    xs = (x[1] + xt)*K.Ry_eV
    i = np ÷ 2
    while E[i] >= xs
        i -= 1
    end
    i -= 1
    while !(E[i] < xs && E[i + 1] >= xs)
        i += 1
    end
    ipos = i
    interp(ip, e) = F[ip] + (F[ip + 1] - F[ip])/(E[ip + 1] - E[ip])*(e - E[ip])
    frate = interp(ipos, xs)*hgh[1]
    for k in 2:m
        e = (x[k] + xt)*K.Ry_eV
        ip = ipos
        while E[ip] < e
            ip += 1
        end
        ip -= 1
        if ip <= np
            ipos = ip
            frate += interp(ipos, e)*hgh[k]
        end
    end
    frate *= delta_rate_coeff
    (; init=idest1, final=nlev, frate=frate, irate=alpha*Float64(lo.g)/Float64(cont.g))
end

# Radiative transition rates from superlevels to spectroscopic levels: r1−rjnd
# = (ne( j), j= 1, jnd); rjnd+1−rjnd+nt = (Te( j), j= 1, jnt);
# rjnd+nt+1−rjnd+nt+nt*nd = ((A( j, j′), j′ = 1, j′nd), j= 1, jnt);
# rnd+nt+nt*nd+1 = λ(Å); i1= nd; i2= nt; i3= i (lower level);
# i4= k (upper level); i5= Z; i6= ionN

# XSTAR data type: 71

const RadiativeSuperDesc = "Transition rates from superlevel to spect. lvls"

struct RadiativeSuper{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    λ::R                # wavelength (Å)
    ne_grid::Vector{R}  # log₁₀ densities (cm⁻³)
    T_grid::Vector{R}   # log₁₀ temperatures (K)
    A::Matrix{R}        # log₁₀ Einstein A (s⁻¹; a single entry above 30 is A itself), size (T, ne)
end

function RadiativeSuper(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Nd, Nt = ivec[1:2]
    RadiativeSuper(Int8(rate), label, Transition(ivec[3], ivec[4]), ivec[5:6]..., rvec[Nd+Nt+Nd*Nt+1],
        Vector(rvec[1:Nd]), Vector(rvec[Nd+1:Nd+Nt]),
        reshape(rvec[Nd+Nt+1:Nd+Nt+Nd*Nt], Int(Nt), Int(Nd)))   # ucalc reads the temperature fastest
end

const superdecay_cap = 1e10               # decay rate cap for ions 96 and 97 (calcium I and II)
const superdecay_capped_ions = (96, 97)
const superdecay_T_margin = 1.0           # log₁₀ T may exceed the table by this
const superdecay_A_linear = 30.0          # a one-point table above this holds A itself
const superdecay_erg_per_eV = Float64(1.602197f-12)   # ucalc's single-precision erg per eV, for the energy

# calt71: the Einstein A of the record at the conditions of the cell and its wavelength
function superlevel_decay(coef::RadiativeSuper, T, den)
    logne, logT = Float64.(coef.ne_grid), Float64.(coef.T_grid)
    logA = Float64.(coef.A)
    nden, ntem = length(logne), length(logT)
    if nden == 1 && ntem == 1
        a = logA[1]
        return (10^(a > superdecay_A_linear ? log10(a) : a), Float64(coef.λ))
    end
    rne = min(log10(den), logne[nden])
    rte = clamp(log10(T), logT[1] - superdecay_T_margin, logT[ntem] + superdecay_T_margin)
    in = 1
    if rne > logne[1]
        in = 0
        while true
            in += 1
            in < nden && rne >= logne[in + 1] && continue
            break
        end
    end
    it = 1
    if rte >= logT[1]
        it = 0
        while true
            it += 1
            if it >= ntem
                it = ntem - 1
            elseif rte >= logT[it + 1]
                continue
            end
            break
        end
    end
    at(j) = logA[it, j] + (logA[it + 1, j] - logA[it, j])/(logT[it + 1] - logT[it])*(rte - logT[it])
    rec = at(in)
    in < nden && (rec += (at(in + 1) - rec)/(logne[in + 1] - logne[in])*(rne - logne[in]))   # (at the last density rne equals logne[in])
    (exp10(rec), Float64(coef.λ))
end

"""
    rate(coef::RadiativeSuper, cell; levels, mass, vturb=1.0, ptmp=(0.5, 0.5), nlev, index=false)

Radiative decay between a superlevel and a spectroscopic level (XSTAR ucalc type 71). The Einstein
A is interpolated in log₁₀ T and log₁₀ n from the table of the record (the temperature may exceed the
table by one dex) and multiplied by the sum `ptmp` of the escape probabilities, with a cap of 10¹⁰ s⁻¹
for ions 96 and 97. As in `ucalc`, the weights come from the first level of the record (`transition.lower`,
`init`) and the decay is returned as `irate` (`frate` is 0). `ienergy` is the decay times the energy of the
line (from its wavelength above 0.1 Å, from the levels otherwise) and `opacity` the line-centre
opacity for the atomic `mass` (amu) and the turbulent velocity `vturb` (km/s), both zero for lines
longer than 10⁹ Å. Needs levels `init` and `final` between 1 and `nlev`.
"""
function rate(coef::RadiativeSuper, cell::Cell; levels, mass=1.0, vturb=1.0, ptmp=(0.5, 0.5),
    nlev=typemax(Int), index=false, verbose=false, kw...)

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0., opacity=0.)
    i1, i2 = coef.transition.lower, coef.transition.upper
    index && return (; none..., init=i1, final=i2)
    (i1 <= 0 || i1 > nlev || i2 <= 0 || i2 > nlev) && return none
    a = get(levels, (coef.ion, i1), nothing)
    b = get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return none
    A, elin = superlevel_decay(coef, cell.T*T_unit, cell.ntot)
    decay = A*(ptmp[1] + ptmp[2])
    coef.ion in superdecay_capped_ions && (decay = min(decay, superdecay_cap))
    flin = A_to_f*A*a.g*elin^2/b.g
    vtherm = sqrt((vturb*cm_per_km)^2 + (thermal_speed/sqrt(mass/cell.T))^2)
    sigma = elin > no_wavelength ? 0. : line_xsec*flin*elin*cm_per_Å/vtherm
    ener = abs(Float64(a.E) - Float64(b.E))
    energy = decay*ener*ergsev
    elin > 0.1 && (energy = decay*hc_eVÅ/(elin + tiny)*superdecay_erg_per_eV)
    (; init=i1, final=i2, frate=0., irate=decay, fenergy=0., ienergy=energy, opacity=sigma)
end

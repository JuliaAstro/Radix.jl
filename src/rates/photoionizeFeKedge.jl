# Photoionization cross sections for Fe ions obtained by summation of
# resonances near the K edge [44]: r1 = Zeff; r2 = Eth (Ryd); r3= f; r4 = γ;
# r5 = scaling factor; i1= n; i2= L; i3= 2J; i4= Z; i5= kN−1; i6= ionN−1;
# i7= iN ; i8= ionN

# XSTAR data type: 85

const PhotoionizeFeKedgeDesc = "Iron K Pi xsections, spectator Auger summed"

struct PhotoionizeFeKedge{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I             # principal quantum number
    L::I             # orbital angular momentum
    twoJ::I          # 2J
    Z::I             # atomic number
    parent::Parent{I}           # parent ion and level
    level::I         # level index
    ion::I           # ion index (XSTAR ionN)
    Zeff::R
    E_th::R          # Ry
    f::R
    γ::R
    scale::R         # scaling factor
end

function PhotoionizeFeKedge(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    PhotoionizeFeKedge(Int8(rate), label, ivec[1:4]..., Parent(ivec[6], ivec[5]), ivec[7:8]..., rvec...)
end

# constants of pexs (resonance series of the Fe K edge) and of the rate
const pexs_nmax = 30                   # highest principal quantum number of the series
const pexs_a_coeff = 8.06725
const pexs_window = 30.0               # series start this many widths below the first resonance
const fe_ion_offset = 114              # ucalc's charge is the ion index minus this
const fe_threshold_fraction = 0.8      # the integral starts at this fraction of the edge energy

"""
    pexs(nmin, zc, eion, far, gam, scal, E)

Photoabsorption cross section (cm², before the 1e-18) of a series of resonances converging
on the K edge, on the energies `E` (Ry), after Palmeri et al. (XSTAR's `pexs`):
`nmin` the lowest resonance, `zc` the effective charge, `eion` the edge energy (Ry), `far`
the oscillator-strength scale, `gam` the resonance width and `scal` a scale factor.
Returns `nothing` if `nmin ≥ 30`.
"""
function pexs(nmin, zc, eion, far, gam, scal, E)
    pexs_pi = constants().pi_pexs
    nmin >= pexs_nmax && return nothing
    x = zeros(pexs_nmax)
    a = zeros(pexs_nmax)
    nres = pexs_nmax
    for n in nmin:pexs_nmax
        x[n] = -(zc/n)^2
        a[n] = pexs_a_coeff*far*nmin^3/n^3
        n > nmin && x[n] - x[n - 1] > gam/2 && (nres = n)
    end
    xmin = x[nmin] - pexs_window*gam
    xres = x[nres]
    kdim = length(E)
    e = E .- eion
    axs = zeros(kdim)
    jmin, jres, jmax = 1, kdim, kdim
    for i in 2:kdim
        e[i - 1] <= xmin && e[i] > xmin && (jmin = i - 1)
        e[i - 1] <= xres && e[i] > xres && (jres = i - 1)
        e[i - 1] < 0 && e[i] >= 0 && (jmax = i - 1)
    end
    jmin == jmax && (jmax = jmin + 1)
    for ii in nmin:pexs_nmax, jj in jmin:jres
        axs[jj] += a[ii]/pexs_pi*gam/2/((e[jj] - x[ii])^2 + (gam/2)^2) +
            a[ii]/pexs_pi*gam/2/((abs(e[jj]) - x[ii])^2 + (gam/2)^2)
    end
    for ij in jres + 1:jmax
        axs[ij] = axs[jres]
    end
    scal .* axs
end

"""
    rate(coef::PhotoionizeFeKedge, cell; radiation, abund=(0, 0), opacity=nothing)

Photoionization of Fe near the K edge from a series of resonances (XSTAR ucalc type 85),
integrated with `photoionization_integrals_fo` from 0.8 of the edge energy. Returns the
photoionization rate `frate` and the heating `fenergy`, `fenergy2`; the recombination
terms and the opacity are zero as in ucalc (the opacity and emissivity arrays are still
filled). `final` is the level 1. `ucalc` computes the charge of the resonance series as
the ion index minus 114, a leftover of an earlier ion numbering, which is reproduced.
"""
function rate(coef::PhotoionizeFeKedge, cell::Cell; radiation, abund=(0.0, 0.0),
    opacity=nothing, ptmp=nothing, lfast=nothing, levels=nothing, nlev=nothing,
    index=false, verbose=false)
    K = constants()

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.,
        fenergy2=0., ienergy2=0., opacity=0.)
    index && return (; none..., init=Int(coef.level), final=1)
    ett2 = Float64(coef.E_th)*K.Ry_eV_single
    σ = pexs(Int(coef.n), Float64(coef.ion - fe_ion_offset), Float64(coef.E_th),
        Float64(coef.f), Float64(coef.γ), Float64(coef.scale), radiation.E ./ K.Ry_eV_single)
    σ === nothing && return none
    r = photoionization_integrals_fo(radiation, ett2*fe_threshold_fraction, σ .* Mb,
        cell.T, 1.0, cell.nₑ; abund=abund, ntot=cell.ntot, lfast=1, opacity=opacity)
    (; init=Int(coef.level), final=1, frate=r.pirt, irate=0., fenergy=r.piht,
        ienergy=0., fenergy2=r.piht2, ienergy2=0., opacity=0.)
end

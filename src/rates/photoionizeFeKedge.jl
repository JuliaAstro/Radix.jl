# Photoionization cross sections for Fe ions obtained by summation of
# resonances near the K edge [44]: r1 = Zeff; r2 = Eth (Ryd); r3= f; r4 = γ;
# r5 = scaling factor; i1= n; i2= L; i3= 2J; i4= Z; i5= kN−1; i6= ionN−1;
# i7= iN ; i8= ionN

# XSTAR data type: 85

const PhotoionizeFeKedgeDesc = "Iron K Pi xsections, spectator Auger summed"

struct PhotoionizeFeKedge{I, R, L<:Levels} <: AbstractRate
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
    levels::L              # the level data of the database (a Levels)
end

function PhotoionizeFeKedge(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    PhotoionizeFeKedge(Int8(rate), label, ivec[1:4]..., Parent(ivec[6], ivec[5]), ivec[7:8]..., rvec..., levels)
end

# constants of pexs (resonance series of the Fe K edge) and of the rate
const pexs_nmax = 30                   # highest principal quantum number of the series
const pexs_a_coeff = 8.06725
const pexs_window = 30.0               # series start this many widths below the first resonance
const fe_ion_offset = 114              # ucalc's charge is the ion index minus this
const fe_threshold_fraction = 0.8      # the integral starts at this fraction of the edge energy

"""
    pexs(nmin, zc, eion, far, gam, scal, E, divisor=1, out)

Photoabsorption cross section (cm², before the 1e-18) of a series of resonances converging
on the K edge, on the energies `E/divisor` (Ry), after Palmeri et al. (XSTAR's `pexs`):
`nmin` the lowest resonance, `zc` the effective charge, `eion` the edge energy (Ry), `far`
the oscillator-strength scale, `gam` the resonance width and `scal` a scale factor.
Returns `nothing` if `nmin ≥ 30`.
"""
function pexs(nmin, zc, eion, far, gam, scal, E, divisor=1.0, out=zeros(length(E)))
    pexs_pi = constants().pi_pexs
    nmin >= pexs_nmax && return nothing
    xr = scratch(:fek_resonances, pexs_nmax)           # the energies and the strengths of the resonances
    ar = scratch(:fek_strengths, pexs_nmax)
    nres = pexs_nmax
    for n in nmin:pexs_nmax
        xr[n] = -(zc/n)^2
        ar[n] = pexs_a_coeff*far*nmin^3/n^3
        n > nmin && xr[n] - xr[n - 1] > gam/2 && (nres = n)
    end
    xmin = xr[nmin] - pexs_window*gam
    xres = xr[nres]
    kdim = length(E)
    e = scratch(:fek_energies, length(E))
    for i in eachindex(E)
        e[i] = E[i]/divisor - eion
    end
    axs = fill!(out, 0.0)
    jmin, jres, jmax = 1, kdim, kdim
    for i in 2:kdim
        e[i - 1] <= xmin && e[i] > xmin && (jmin = i - 1)
        e[i - 1] <= xres && e[i] > xres && (jres = i - 1)
        e[i - 1] < 0 && e[i] >= 0 && (jmax = i - 1)
    end
    jmin == jmax && (jmax = jmin + 1)
    for jj in jmin:jres                                # (each energy sums the resonances in the same order, in a local sum: the work arrays of the task may alias one another, which would keep the loop from being optimized)
        ej, sum = e[jj], axs[jj]
        for ii in nmin:pexs_nmax
            sum += ar[ii]/pexs_pi*gam/2/((ej - xr[ii])^2 + (gam/2)^2) +
                ar[ii]/pexs_pi*gam/2/((abs(ej) - xr[ii])^2 + (gam/2)^2)
        end
        axs[jj] = sum
    end
    for ij in jres + 1:jmax
        axs[ij] = axs[jres]
    end
    axs .*= scal
end

"""
    rate(coef::PhotoionizeFeKedge, cell; abund=(0, 0), opacity=nothing)

Photoionization of Fe near the K edge from a series of resonances (XSTAR ucalc type 85),
integrated with `photoionization_integrals_fo` from 0.8 of the edge energy. Returns the
photoionization rate `frate` and the heating `fenergy`, `fenergy2`; the recombination
terms and the opacity are zero as in ucalc (the opacity and emissivity arrays are still
filled). `final` is the level 1. `ucalc` computes the charge of the resonance series as
the ion index minus 114, a leftover of an earlier ion numbering, which is reproduced.
"""
function rate(coef::PhotoionizeFeKedge, cell::Cell; radiation=NO_RADIATION, abund=(0.0, 0.0), opacity=nothing, ptmp=nothing, lfast=nothing, index=false)
    K = constants()

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.,
        fenergy2=0., ienergy2=0., opacity=0.)
    index && return (; none..., init=Int(coef.level), final=1)
    ett2 = Float64(coef.E_th)*K.Ry_eV_single
    σ = pexs(Int(coef.n), Float64(coef.ion - fe_ion_offset), Float64(coef.E_th),
        Float64(coef.f), Float64(coef.γ), Float64(coef.scale), radiation.E, K.Ry_eV_single, scratch(:fek_sigma, length(radiation.E)))
    σ === nothing && return none
    σ .*= Mb
    r = photoionization_integrals_fo(radiation, ett2*fe_threshold_fraction, σ,
        cell.T, 1.0, cell.nₑ; abund=abund, ntot=cell.ntot, lfast=1, opacity=opacity)
    (; init=Int(coef.level), final=1, frate=r.pirt, irate=0., fenergy=r.piht,
        ienergy=0., fenergy2=r.piht2, ienergy2=0., opacity=0.)
end

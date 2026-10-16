# Hydrogenic photoionization cross sections, in the format of Bautista: i1= n (principal quantum number); i2= l (orbital quantum number); i3= the charge of the ion; i4.. = not used;
# the last two: the level i and ionN; the reals are not used

# XSTAR data type: 64

const PhotoionizeHydrogenicDesc = "hydrogenic pi xsections, bautista format"

struct PhotoionizeHydrogenic{I, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I        # principal quantum number
    l::I        # orbital angular momentum
    charge::I   # charge of the ion
    level::I    # level index
    ion::I      # ion index (XSTAR ionN)
    levels::L              # the level data of the database (a Levels)
end

function PhotoionizeHydrogenic(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    PhotoionizeHydrogenic(Int8(rate), label, ivec[1], ivec[2], ivec[3], ivec[end - 1], ivec[end], levels)
end

const hydrogenic_ett_min = 1e-5          # eV: the levels with a smaller ionization potential give no rates
const hydrogenic_Ry = Float64(13.605692f0)
const hydrogenic_Ry_milne = Float64(13.6f0)

"""
    rate(coef::PhotoionizeHydrogenic, cell; radiation=NO_RADIATION, abund=(0, 0), opacity=nothing, ptmp=nothing, lfast=nothing, index=false)

Photoionization of a level of the principal quantum number `n` and angular momentum `l` of an ion of charge `charge` with the hydrogenic cross sections of Brocklehurst (1971) (`hphotx`) (XSTAR ucalc type 64), integrated over `radiation`
with the cross section on the grid itself (`photoionization_integrals_fo`) from the threshold, the energy `|E_c − E|` of the level to the continuum level, to the top of the grid: `frate` is the photoionization rate and `fenergy`,
`fenergy2` the heating. `irate` is the recombination coefficient of the Milne relation of the same cross section times the weight of the level (the weight of the continuum level is not used; ucalc does not multiply by
the electron density), and `ienergy`, `ienergy2` the recombination cooling of the integrals. `init` is the level and `final` the continuum level (`nlev`). The records with a threshold below 10⁻⁵ eV give no rates. ucalc
takes the top of the grid from the value that an earlier call left (`nphint`): the last bin that `enxt` gives, `ncn2 - max(2, ncn2/50)`, is taken here.
"""
function rate(coef::PhotoionizeHydrogenic, cell::Cell; radiation=NO_RADIATION, abund=(0.0, 0.0), opacity=nothing, ptmp=nothing, lfast=nothing, index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    none = (; init=Int32(0), final=Int32(0), frate=0., irate=0., fenergy=0., ienergy=0., fenergy2=0., ienergy2=0., opacity=0.)
    idest1 = Int(coef.level)
    index && return (; none..., init=idest1, final=nlev)
    lo = get(levels, (coef.ion, idest1), nothing)
    cont = get(levels, (coef.ion, nlev), nothing)
    (lo === nothing || cont === nothing) && return none
    ett = abs(Float64(cont.E) - Float64(lo.E))
    ett <= hydrogenic_ett_min && return none
    ncn2 = length(radiation.E)
    nphint = ncn2 - max(grid_min_bins, ncn2 ÷ grid_guard_fraction)
    nb1 = nbin(radiation, ett)
    swrat = Float64(lo.g)
    sg = fill!(scratch(:fo_cross_section, ncn2), 0.0)
    erel = scratch(:hydrogenic_energies, max(nphint - nb1 + 1, 0))
    xs = scratch(:hydrogenic_sections, max(nphint - nb1 + 1, 0))
    gu, gl = scratch(:gull_u, max(Int(coef.n), 3) + 1), scratch(:gull_l, max(Int(coef.n), 3) + 1)
    mm = 0
    for ll in nb1:nphint
        mm += 1
        e = max(0.0, (radiation.E[ll] - ett)/hydrogenic_Ry)
        x = hphotx!(gu, gl, e, Int(coef.charge), Int(coef.n), Int(coef.l))
        sg[ll] = x*Mb
        xs[mm] = x
        erel[mm] = e
    end
    r = photoionization_integrals_fo(radiation, ett, sg, cell.T, swrat, cell.nₑ; abund=abund, ntot=cell.ntot, lfast=1, opacity=opacity)
    α = milne_recombination(cell.T*T_unit, view(erel, 1:mm), view(xs, 1:mm), ett/hydrogenic_Ry_milne)*swrat
    (; init=idest1, final=nlev, frate=r.pirt, irate=α, fenergy=r.piht, ienergy=r.rrcl, fenergy2=r.piht2, ienergy2=r.rrcl2, opacity=r.opakab)
end

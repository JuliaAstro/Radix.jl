# Partial photoionization cross section of iN -th level of the N-electron ion
# leaving the (N−1)-electron ion in the kN−1-th level [51]: r1 = E(th) (eV);
# r2 = E(0) (eV); r3= σ(0) (Mb); r4 = y(a); r5 = P; r6= y(w) ; i1= N;
# i2= n (shell principal quantum number); i3= l (orbital quantum number of the
# subshell); i4= kN−1; i5= ionN−1; i6= iN ; i7= ionN ; s1= shell-ion identifier

# XSTAR data type: 59

const ParPhotoIonize3Desc = "verner pi x!"

@with_levels struct ParPhotoIonize3{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n_electrons::I   # electrons in the ion
    n::I             # shell principal quantum number
    l::I             # subshell orbital quantum number
    parent::Parent{I}           # parent ion and level
    level::I         # level index
    ion::I           # ion index (XSTAR ionN)
    E_th::R          # eV
    E0::R            # eV
    σ0::R            # Mb
    ya::R
    P::R
    yw::R
end

function ParPhotoIonize3(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ParPhotoIonize3(Int8(rate), label, ivec[1:3]..., Parent(ivec[5], ivec[4]), ivec[6:7]..., rvec...)
end

const ph3_qq_base = 5.5               # exponent of the fit: 5.5 + l - P/2
const ph3_exp_limit = 60.0
const ph3_floor = 1e-48
const ph3_flux_min = 1e-20            # no rate if the integrated spectrum above the threshold is below this

"""
    rate(coef::ParPhotoIonize3, cell; abund=(0, 0), opacity=nothing)

Photoionization from the analytic fit of Verner et al. (XSTAR ucalc type 59),
evaluated on the energy grid of `radiation` and integrated with `photoionization_integrals_fo`
(the recombination terms are zero for rate type 1 and for excited levels, as in ucalc).
The levels and their number come from the coefficient's level table (`attach_levels`). Returns `init` and
`final` (`nlev` + parent level − 1), `frate` and `irate`, the energies `fenergy`, `ienergy`,
`fenergy2`, `ienergy2` and the `opacity`. Nothing is returned if the parent ion is more
than one step away, the threshold is outside the grid, or the integrated spectrum above it
is negligible. Only the 6-real form of the record is supported. `ptmp` and `lfast` are
accepted for uniformity with the other photoionization rates but not used.
"""
function rate(coef::ParPhotoIonize3, cell::Cell; radiation=NO_RADIATION, abund=(0.0, 0.0), opacity=nothing, ptmp=nothing, lfast=nothing, index=false)
    levels = levels_of(coef)
    nlev = nlevels(levels, coef.ion)

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.,
        fenergy2=0., ienergy2=0., opacity=0.)
    coef.parent.ion > coef.ion + 1 && return none
    idest1 = Int(coef.level)
    idest2 = max(nlev + Int(coef.parent.level) - 1, 1)
    index && return (; none..., init=idest1, final=idest2)
    (idest1 > nlev || idest1 <= 0) && return none
    cont = get(levels, (coef.ion, nlev), nothing)
    ground = get(levels, (coef.ion, 1), nothing)
    (cont === nothing || ground === nothing) && return none
    ncn2 = length(radiation.E)
    numcon2 = max(grid_min_bins, ncn2 ÷ grid_guard_fraction)
    ett_ion = Float64(cont.E)                       # ionization energy of the ion
    ett_ion <= 0 && return none
    nphint = max(ncn2 - numcon2, nbin(radiation, ett_ion) + 1)
    nbin(radiation, ett_ion) >= nphint - 1 && return none
    ggup = Float64(cont.g)
    if idest2 > nlev
        par = get(levels, (coef.parent.ion, Int(coef.parent.level)), nothing)
        par === nothing && return none
        ggup = Float64(par.g)
    end
    ggup <= min_g && return none
    swrat = Float64(ground.g)/ggup
    ett = Float64(coef.E_th)
    nb1 = nbin(radiation, ett)
    nb1 >= nphint - 1 && return none
    radiation.Fint[nb1] < ph3_flux_min && return none

    E0, s0, ya, pp, yw = Float64(coef.E0), Float64(coef.σ0), Float64(coef.ya),
        Float64(coef.P), Float64(coef.yw)
    qq = ph3_qq_base + coef.l - pp/2
    sg = zeros(ncn2)
    for ll in nb1:ncn2 - numcon2
        xx = radiation.E[ll]/E0
        yyqq = exp(-min(ph3_exp_limit, max(-ph3_exp_limit, qq*log(max(ph3_floor, xx)))))
        sg[ll] = s0*((xx - 1)^2 + yw^2)*yyqq*(1 + sqrt(xx/ya))^(-pp)*Mb
    end
    r = photoionization_integrals_fo(radiation, ett, sg, cell.T, swrat, cell.nₑ;
        abund=abund, ntot=cell.ntot, lfast=1, opacity=opacity)
    recomb = !(coef.rtype == 1 || idest1 > 1)
    (; init=idest1, final=idest2, frate=r.pirt, irate=recomb ? r.rrrt : 0.,
        fenergy=r.piht, ienergy=recomb ? r.rrcl : 0., fenergy2=r.piht2,
        ienergy2=recomb ? r.rrcl2 : 0., opacity=r.opakab)
end

# Integrals of photoionization cross sections over the radiation field:
# XSTAR's phextrap() and phint53() (xstarlib/src).

const phextrap_step = 1.3            # energy factor between extrapolated points
const phextrap_floor = 1e-27         # stop extrapolating below this cross section (cm²)
const phextrap_Emax = 2e5            # ... or above this energy (eV)
const phextrap_Ry = 13.6             # eV per Ry in the extrapolation

"""
    phextrap(ε, σ, E_th, ncn2)

Extends a cross-section table (`ε` in Ry above the threshold `E_th` eV, `σ` in
cm²) with points falling as E⁻³ up to 2e5 eV, as XSTAR does. The last tabulated
point is replaced by the first extrapolated one (or dropped if there is none).
"""
function phextrap(ε::AbstractVector, σ::AbstractVector, E_th, ncn2)
    n = length(ε)
    n <= 0 && return (Float64.(ε), Float64.(σ))
    eout = Float64.(ε[1:n-1]); sout = Float64.(σ[1:n-1])
    n < 2 && return (eout, sout)
    s1 = σ[n-1]
    e1 = ε[n-1]*phextrap_Ry + E_th
    nadd = 0
    while s1 > phextrap_floor && nadd + n < ncn2 && e1 < phextrap_Emax
        e2 = e1*phextrap_step
        s2 = s1/phextrap_step^3
        nadd += 1
        push!(sout, s2)
        push!(eout, (e2 - E_th)/phextrap_Ry)
        e1, s1 = e2, s2
    end
    (eout, sout)
end

# constants of phint53
const ph_min_dE = 1e-8               # eV; flat segments below this width
const ph_tiny_dE = 1e-36
const ph_tiny = 1e-24
const bb_energy_cap = 2e4            # eV; the Planck-like factor stops growing here
const bb_coeff = 1.571e22
const ph_exp_limit = 200.0           # recombination terms only while (E - E_th)/kT is below this

"""
    photoionization_integrals(rad, E_th, ε, σ, T, rnist, ptmp; abund=(0, 0), ntot=0, lfast=1, opacity=nothing)

Photoionization rate `pirt`, recombination rate `rrrt`, and the energy-weighted
`piht`, `rrcl` (erg s⁻¹, measured from zero) and `piht2`, `rrcl2` (measured from
the threshold), by integrating the cross section table (`ε` in Ry above the
threshold `E_th` in eV, `σ` in cm²) over the radiation field (XSTAR's `phint53`).
`T` is in 10⁴ K, `rnist` the Saha factor and `ptmp` the two escape
probabilities (in the reverse and forward direction). The recombination terms need `lfast ≥ 2`. If `opacity` is given,
the continuum opacity and recombination emissivity are added to it. Also returns
`opakab`, the opacity of the first bins above the threshold.
"""
function photoionization_integrals(rad::Radiation, E_th, ε, σ, T, rnist, ptmp;
    abund=(0.0, 0.0), ntot=0.0, lfast=1, opacity=nothing)

    epi, bremsa = rad.E, rad.F
    ncn2 = length(epi)
    ntmp = length(ε)
    out = (; pirt=0.0, rrrt=0.0, piht=0.0, rrcl=0.0, piht2=0.0, rrcl2=0.0, opakab=0.0)
    ntmp <= 0 && return out
    abund1, abund2 = abund
    pesc = ptmp[1] + ptmp[2]
    eth = E_th
    bktm = kB_cgs*T*T_unit/ergsev
    nphint = ncn2 - max(grid_min_bins, ncn2 ÷ grid_guard_fraction)
    sumr = sumh = sumh2 = sumc = sumc2 = sumi = 0.0
    ener = eth + ε[1]*Ry_eV
    nb1 = nbin(rad, ener)
    while epi[nb1] < ener && nb1 < nphint
        nb1 += 1
    end
    nb1 = max(nb1 - 1, 1)
    nb1 >= nphint && return out

    enermx = eth + ε[ntmp]*Ry_eV
    nbn = max(nbin(rad, enermx), min(nb1 + 1, ncn2 - 1))
    sgbar = zeros(ncn2 + 1)
    kl = nb1
    jk = 1
    e1 = epi[kl]
    e2 = eth + ε[jk]*Ry_eV
    s2 = σ[jk]
    if e1 < e2
        kl += 1
        e1 = epi[kl]
    end
    e1o = e2
    s2o, e2o, s2t = s2, e2, s2        # XSTAR leaves these undefined if the table has < 3 points
    sum = 0.0
    done = false
    while !done
        while e2 < e1 && jk < ntmp - 1
            jk += 1
            e2o = e2
            s2o = s2
            e2 = eth + ε[jk]*Ry_eV
            s2 = σ[jk]
            sum += (s2 + s2o)*(e2 - e2o)/2
        end
        sum -= (s2 + s2o)*(e2 - e2o)/2
        e2t = e1
        s2t = e2 - e2o > ph_min_dE ? s2o + (s2 - s2o)*(e2t - e2o)/(e2 - e2o + ph_tiny) : s2o
        sum += (s2t + s2o)*(e2t - e2o)/2
        sgbar[kl] = abs(e1 - e1o) > ph_tiny_dE ? sum/(e1 - e1o) : 0.0
        e1o = e1
        kl += 1
        e1 = epi[kl]
        while e1 < e2 && kl < ncn2
            e2t = e1
            s2t = e2 - e2o > ph_min_dE ? s2o + (s2 - s2o)*(e2t - e2o)/(e2 - e2o) : s2o
            sum = (s2t + s2t)*(e1 - e1o)/2
            sgbar[kl] = sum/(e1 - e1o)
            e1o = e1
            kl += 1
            e1 = epi[kl]
        end
        sum = (s2 + s2t)*(e2 - e2t)/2
        (kl > nphint - 1 || jk >= ntmp - 1) && (done = true)
    end
    klmax = kl - 1

    bb(e) = min(bb_energy_cap, e)^3*bb_coeff*2
    sgtpp = sgbar[nb1]
    bremtmpp = bremsa[nb1]/fourpi_xstar
    epiip = epi[nb1]
    temprp = fourpi_xstar*sgtpp*bremtmpp/epiip
    temphp = temprp*epiip
    temphp2 = temprp*(epiip - eth)
    exptst = (epiip - eth)/bktm
    exptmpp = expo(-exptst)
    tempip = rnist*bb(epiip)*sgtpp*exptmpp/epiip*pesc
    tempcp = tempip*epiip
    tempcp2 = tempip*(epiip - eth)
    tempi = tempc = tempc2 = 0.0
    opakab = 0.0
    kl = nb1
    while kl < klmax
        sgtp = max(0.0, sgbar[kl])
        sgtpp = sgbar[kl + 1]
        bremtmpp = bremsa[kl + 1]/fourpi_xstar
        epii = epi[kl]
        epiip = epi[kl + 1]
        tempr = temprp
        temprp = fourpi_xstar*sgtpp*bremtmpp/epiip
        wwir = (epiip - epii)/2
        sumr += tempr*wwir + temprp*wwir
        temph, temph2 = temphp, temphp2
        temphp = temprp*epiip
        temphp2 = temprp*(epiip - eth)
        sumh += temph*wwir + temphp*wwir
        sumh2 += temph2*wwir + temphp2*wwir
        exptsto = exptst
        exptst = (epiip - eth)/bktm
        exptmpp = 0.0
        if exptsto < ph_exp_limit && lfast >= 2
            exptmpp = expo(-exptst)
            tempi = tempip
            tempip = rnist*bb(epiip)*sgtpp*exptmpp*fourpi_xstar/epiip
            atmp2 = tempip*epiip
            tempip *= pesc
            sumi += tempi*wwir + tempip*wwir
            tempc, tempc2 = tempcp, tempcp2
            tempcp = tempip*epiip
            tempcp2 = tempip*(epiip - eth)
            sumc += tempc*wwir + tempcp*wwir
            sumc2 += tempc2*wwir + tempcp2*wwir
            if opacity !== nothing
                opacity.emissivity[1, kl] += abund2*atmp2*ptmp[1]*ntot/fourpi_xstar
                opacity.emissivity[2, kl] += abund2*atmp2*ptmp[2]*ntot/fourpi_xstar
            end
        end
        optmp = abund1*ntot*sgtp
        if opacity !== nothing
            opacity.total[kl] += optmp
            opacity.continuum[kl] += optmp
        end
        if kl == nb1 + 2
            optmp2 = rnist*exptmpp*sgtp*abund2*pesc*ntot
            opakab = max(0.0, optmp - optmp2)
        end
        kl += 1
    end
    (; pirt=sumr, rrrt=sumi, piht=sumh*ergsev, rrcl=sumc*ergsev,
        piht2=sumh2*ergsev, rrcl2=sumc2*ergsev, opakab=opakab)
end

const saha_coeff = 2.07e-16            # Saha-Boltzmann factor (cm³ K^{3/2})
const saha_T_exp = -1.5
const min_g = 1e-24                    # statistical weights below this are ignored
const rnisseu_floor = 1e-37
const heat_floor = 1e-43


"""
    photoionize_level(coef, cell; levels, radiation, nlev, ptmp=(0.5, 0.5), abund=(0, 0), lfast=1, opacity=nothing, extrapolate, shifted, parent=true, rates_only=false)

Photoionization of a level of an ion by integrating the cross-section table of
`coef` (fields `level`, `ion`, `parent`, `E_grid` in Ry above the threshold and `σ`
in Mb) over `radiation` (XSTAR ucalc types 49 and 53). `levels` is a
`level_table`, `nlev` the number of levels of the ion (the last is the
continuum, the ground state of the next ion), `ptmp` the two escape
probabilities and `abund` the populations of the initial and final levels (for the
opacity). `lfast ≥ 2` also integrates the recombination terms. If `opacity` is an
`Opacity` the continuum opacity and recombination emissivity are added to it.
`extrapolate` extends the table with an E⁻³ tail (types 49 and 88) and `shifted` adds the
excitation energy of an excited parent level to the threshold (type 53). With
`parent=false` the final level is always the continuum (type 88), and
`rates_only=true` keeps only the photoionization rate and the opacity, zeroing the
rest as ucalc type 88 does (the opacity and emissivity arrays are still filled).

Returns `init` (the level) and `final` (`nlev` + parent level − 1), `frate` the
photoionization rate and `irate` the radiative recombination rate (s⁻¹),
`fenergy`/`ienergy` the heating and recombination cooling (erg s⁻¹), `fenergy2`
and `ienergy2` the same measured from the energies of the two levels, and
`opacity` the opacity of the first bins above the threshold. For parent levels
above the ground state, `ucalc` reads the energy of the final level from stale
memory; here it is the continuum energy plus the energy of the parent level.
"""
function photoionize_level(coef, cell::Cell; levels, radiation, nlev,
    ptmp=(0.5, 0.5), abund=(0.0, 0.0), lfast=1, opacity=nothing, index=false,
    extrapolate, shifted, parent=true, rates_only=false)

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.,
        fenergy2=0., ienergy2=0., opacity=0.)
    idest1 = Int(coef.level)
    idest2 = parent ? nlev + max(0, Int(coef.parent.level)) - 1 : nlev
    (idest1 >= nlev || idest1 <= 0) && return none
    isempty(coef.E_grid) && return none                  # no cross-section table
    index && return (; none..., init=idest1, final=idest2)
    lo = get(levels, (coef.ion, idest1), nothing)
    cont = get(levels, (coef.ion, nlev), nothing)
    (lo === nothing || cont === nothing) && return none
    eth = Float64(lo.E_inf) - Float64(lo.E)
    eth <= 0 && return none

    ε = Float64.(coef.E_grid)
    σ = max.(Float64.(coef.σ)*Mb, 0.0)
    ggup = Float64(cont.g)
    e2 = Float64(cont.E)
    if idest2 > nlev
        par = get(levels, (coef.parent.ion, Int(coef.parent.level)), nothing)
        par === nothing && return none
        ggup = Float64(par.g)
        e2 += Float64(par.E)
        shifted && (eth += Float64(par.E))
    end
    ggup <= min_g && return none
    ε0 = ε[1]
    if extrapolate
        ε, σ = phextrap(ε, σ, eth, length(radiation.E))
        isempty(ε) && return none
    end

    T = cell.T
    q2 = saha_coeff*cell.nₑ*(T*T_unit)^saha_T_exp
    rnissel = Float64(lo.g)*q2/Float64(cont.g)
    ethtmp = shifted ? max(0.0, eth - Float64(cont.E)) : 0.0
    rnist = rnissel*exp(-max(0.0, ethtmp + Ry_eV*ε0)/kT_eV/T)/(rnisseu_floor + 1)
    r = photoionization_integrals(radiation, eth, ε, σ, T, rnist, ptmp;
        abund=abund, ntot=cell.ntot, lfast=lfast, opacity=opacity)

    rates_only && return (; none..., init=idest1, final=idest2, frate=r.pirt,
        opacity=r.opakab)

    dE = abs(e2 - Float64(lo.E))
    fenergy2 = r.piht2*(r.piht - dE*ergsev*r.pirt)/max(heat_floor, r.piht - eth*ergsev*r.pirt)
    ienergy2 = r.rrcl2*(r.rrcl - dE*ergsev*r.rrrt)/max(heat_floor, r.rrcl - eth*ergsev*r.rrrt)
    (; init=idest1, final=idest2, frate=r.pirt, irate=r.rrrt, fenergy=r.piht,
        ienergy=r.rrcl, fenergy2=fenergy2, ienergy2=ienergy2, opacity=r.opakab)
end

# ---------------------------------------------------------------------------
# phintfo and enxt: integration with a cross section given on the grid itself

const enxt_eth_factor = 3.0          # lfast = 3: integrate up to max(3 E_th, E_th + 3 kT) ...
const enxt_kT_factor = 3.0
const enxt_samples = 16              # ... in about this many steps
const enxt_Emax = 1e4                # lfast ≥ 4: integrate up to this energy (eV)

"""
    next_bin(rad, E_th, nb1, T, lfast)

Step `nskp` to the next sampled bin of the photoionization integrals and the last bin
`nphint` to use (XSTAR's `enxt`).
"""
function next_bin(rad::Radiation, E_th, nb1, T, lfast)
    ncn2 = length(rad.E)
    bktm = kB_cgs*T*T_unit/ergsev
    numcon2 = max(grid_min_bins, ncn2 ÷ grid_guard_fraction)
    if lfast <= 2
        nphint = ncn2 - numcon2
        nskp = 1
    elseif lfast == 3
        nphint = nbin(rad, max(enxt_eth_factor*E_th, E_th + enxt_kT_factor*bktm))
        nphint = max(nphint, nb1 + 1)
        nskp = max(1, (nphint - nb1) ÷ enxt_samples)
    else
        nphint = nbin(rad, enxt_Emax)
        nskp = 1
    end
    nphint = max(nphint, nb1 + nskp)
    (nskp, min(nphint, ncn2 - numcon2))
end

const fo_saha = 5.216e-21            # Saha-Boltzmann factor of phintfo (cm³ K^{3/2})
const fo_fourpi = 25.3               # phintfo's value of 8π
const fo_exptst_floor = 1e-36

"""
    photoionization_integrals_fo(rad, E_th, σ, T, swrat, xnx; abund=(0, 0), ntot=0, lfast=1, opacity=nothing)

Photoionization and recombination integrals for a cross section `σ` (cm²) given on the
energy grid of `rad` (XSTAR's `phintfo`): `pirt`, `rrrt`, `piht`, `rrcl`, `piht2`,
`rrcl2` as in `photoionization_integrals`, and the opacity `opakab` of the first bins.
`swrat` is the ratio of the statistical weights of the two levels and `xnx` the
electron density. The recombination terms always run (the convergence test of the
original never stops the loop).
"""
function photoionization_integrals_fo(rad::Radiation, E_th, σ, T, swrat, xnx;
    abund=(0.0, 0.0), ntot=0.0, lfast=1, opacity=nothing)

    epi, bremsa = rad.E, rad.F
    ncn2 = length(epi)
    abund1 = abund[1]
    eth = E_th
    nb1 = nbin(rad, eth)
    bktm = kB_cgs*T*T_unit/ergsev
    rnist = fo_saha*swrat/T/sqrt(T)
    sumr = sumh = sumh2 = sumc = sumc2 = sumi = 0.0
    tempr = tempi = atmp2 = atmp22 = 0.0
    opakab = 0.0
    nphint = max(ncn2 - max(grid_min_bins, ncn2 ÷ grid_guard_fraction), nb1 + 1)
    ener = epi[nb1]
    kl = nb1
    while kl <= nphint
        enero = ener
        ener = epi[kl]
        epii = ener
        sgtmp = σ[kl]
        bremtmp = bremsa[kl]/fo_fourpi
        tempro = tempr
        tempr = fo_fourpi*sgtmp*bremtmp/epii
        deld = ener - enero
        sumr += (tempr + tempro)*deld/2
        sumh += (tempr*ener + tempro*enero)*deld*ergsev/2
        sumh2 += (tempr*(ener - eth) + tempro*(enero - eth))*deld*ergsev/2
        exptmp = expo(-max(fo_exptst_floor, (epii - eth)/bktm))
        bbnurj = min(epii, bb_energy_cap)^3*bb_coeff
        tempi1 = rnist*bbnurj*exptmp*sgtmp/epii
        tempi2 = rnist*bremtmp*exptmp*sgtmp/epii
        tempio = tempi
        tempi = tempi1 + tempi2
        atmp2o = atmp2
        atmp2 = tempi1*epii
        atmp22o = atmp22
        atmp22 = tempi1*(epii - eth)
        sumi += (tempi + tempio)*deld/2
        sumc += (atmp2 + atmp2o)*deld*ergsev/2
        sumc2 += (atmp22 + atmp22o)*deld*ergsev/2
        optmp = abund1*sgtmp*ntot
        kl <= nb1 + 1 && (opakab = optmp)
        if opacity !== nothing
            opacity.total[kl] += optmp
            opacity.continuum[kl] += optmp
        end
        nskp, nphint = next_bin(rad, eth, nb1, T, lfast)
        kl += nskp
    end
    (; pirt=sumr, rrrt=xnx*sumi, piht=sumh, rrcl=xnx*sumc, piht2=sumh2,
        rrcl2=xnx*sumc2, opakab=opakab)
end

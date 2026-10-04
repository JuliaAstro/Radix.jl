# Integrals of photoionization cross sections over the radiation field:
# XSTAR's phextrap() and phint53() (xstarlib/src).

const phextrap_step = 1.3            # energy factor between extrapolated points
const phextrap_floor = 1e-27         # stop extrapolating below this cross section (cm²)
const phextrap_Emax = 2e5            # ... or above this energy (eV)

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
    Ry = constants().Ry_eV_coarse
    e1 = ε[n-1]*Ry + E_th
    nadd = 0
    while s1 > phextrap_floor && nadd + n < ncn2 && e1 < phextrap_Emax
        e2 = e1*phextrap_step
        s2 = s1/phextrap_step^3
        nadd += 1
        push!(sout, s2)
        push!(eout, (e2 - E_th)/Ry)
        e1, s1 = e2, s2
    end
    (eout, sout)
end

# constants of phint53
const ph_min_dE = 1e-8               # eV; flat segments below this width
const ph_tiny_dE = 1e-36
const ph_tiny = 1e-24
const bb_energy_cap = 2e4            # eV; the Planck-like factor stops growing here
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
    K = constants()

    epi, bremsa = rad.E, rad.F
    ncn2 = length(epi)
    ntmp = length(ε)
    out = (; pirt=0.0, rrrt=0.0, piht=0.0, rrcl=0.0, piht2=0.0, rrcl2=0.0, opakab=0.0)
    ntmp <= 0 && return out
    abund1, abund2 = abund
    pesc = ptmp[1] + ptmp[2]
    eth = E_th
    bktm = K.kB_cgs*T*T_unit/K.ergsev
    nphint = ncn2 - max(grid_min_bins, ncn2 ÷ grid_guard_fraction)
    sumr = sumh = sumh2 = sumc = sumc2 = sumi = 0.0
    ener = eth + ε[1]*K.Ry_eV
    nb1 = nbin(rad, ener)
    while epi[nb1] < ener && nb1 < nphint
        nb1 += 1
    end
    nb1 = max(nb1 - 1, 1)
    nb1 >= nphint && return out

    enermx = eth + ε[ntmp]*K.Ry_eV
    nbn = max(nbin(rad, enermx), min(nb1 + 1, ncn2 - 1))
    sgbar = zeros(promote_type(eltype(ε), eltype(σ), typeof(eth)), ncn2 + 1)
    kl = nb1
    jk = 1
    e1 = epi[kl]
    e2 = eth + ε[jk]*K.Ry_eV
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
            e2 = eth + ε[jk]*K.Ry_eV
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

    bb(e) = min(bb_energy_cap, e)^3*K.bb_coeff*2
    sgtpp = sgbar[nb1]
    bremtmpp = bremsa[nb1]/K.fourpi
    epiip = epi[nb1]
    temprp = K.fourpi*sgtpp*bremtmpp/epiip
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
        bremtmpp = bremsa[kl + 1]/K.fourpi
        epii = epi[kl]
        epiip = epi[kl + 1]
        tempr = temprp
        temprp = K.fourpi*sgtpp*bremtmpp/epiip
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
            tempip = rnist*bb(epiip)*sgtpp*exptmpp*K.fourpi/epiip
            atmp2 = tempip*epiip
            tempip *= pesc
            sumi += tempi*wwir + tempip*wwir
            tempc, tempc2 = tempcp, tempcp2
            tempcp = tempip*epiip
            tempcp2 = tempip*(epiip - eth)
            sumc += tempc*wwir + tempcp*wwir
            sumc2 += tempc2*wwir + tempcp2*wwir
            if opacity !== nothing
                opacity.emissivity[1, kl] += abund2*atmp2*ptmp[1]*ntot/K.fourpi
                opacity.emissivity[2, kl] += abund2*atmp2*ptmp[2]*ntot/K.fourpi
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
    (; pirt=sumr, rrrt=sumi, piht=sumh*K.ergsev, rrcl=sumc*K.ergsev,
        piht2=sumh2*K.ergsev, rrcl2=sumc2*K.ergsev, opakab=opakab)
end

const saha_T_exp = -1.5
const min_g = 1e-24                    # statistical weights below this are ignored
const rnisseu_floor = 1e-37
const heat_floor = 1e-43


"""
    photoionize_level(coef, cell; radiation=NO_RADIATION, ptmp=(0.5, 0.5), abund=(0, 0), lfast=1, opacity=nothing, extrapolate=false, shifted=false, parent=true, rates_only=false)

Photoionization of a level of an ion by integrating the cross-section table of
`coef` (fields `level`, `ion`, `parent`, `E_grid` in Ry above the threshold and `σ`
in Mb) over `radiation` (XSTAR ucalc types 49 and 53). The level data, with the number of levels of
the ion (the last is the continuum, the ground state of the next ion), comes from the coefficient
(`attach_levels`); `ptmp` is the two escape probabilities and `abund` the populations of the initial and final levels (for the
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
function photoionize_level(coef, cell::Cell; radiation=NO_RADIATION,
    ptmp=(0.5, 0.5), abund=(0.0, 0.0), lfast=1, opacity=nothing, index=false,
    extrapolate=false, shifted=false, parent=true, rates_only=false)
    K = constants()
    levels = levels_of(coef)
    nlev = nlevels(levels, coef.ion)

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
    q2 = K.saha_coeff*cell.nₑ*(T*T_unit)^saha_T_exp
    rnissel = Float64(lo.g)*q2/Float64(cont.g)
    ethtmp = shifted ? max(0.0, eth - Float64(cont.E)) : 0.0
    rnist = rnissel*exp(-max(0.0, ethtmp + K.Ry_eV*ε0)/K.kT_eV/T)/(rnisseu_floor + 1)
    r = photoionization_integrals(radiation, eth, ε, σ, T, rnist, ptmp;
        abund=abund, ntot=cell.ntot, lfast=lfast, opacity=opacity)

    rates_only && return (; none..., init=idest1, final=idest2, frate=r.pirt,
        opacity=r.opakab)

    dE = abs(e2 - Float64(lo.E))
    fenergy2 = r.piht2*(r.piht - dE*K.ergsev*r.pirt)/max(heat_floor, r.piht - eth*K.ergsev*r.pirt)
    ienergy2 = r.rrcl2*(r.rrcl - dE*K.ergsev*r.rrrt)/max(heat_floor, r.rrcl - eth*K.ergsev*r.rrrt)
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
    K = constants()
    ncn2 = length(rad.E)
    bktm = K.kB_cgs*T*T_unit/K.ergsev
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
    K = constants()

    epi, bremsa = rad.E, rad.F
    ncn2 = length(epi)
    abund1 = abund[1]
    eth = E_th
    nb1 = nbin(rad, eth)
    bktm = K.kB_cgs*T*T_unit/K.ergsev
    rnist = K.fo_saha*swrat/T/sqrt(T)
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
        bremtmp = bremsa[kl]/K.eightpi
        tempro = tempr
        tempr = K.eightpi*sgtmp*bremtmp/epii
        deld = ener - enero
        sumr += (tempr + tempro)*deld/2
        sumh += (tempr*ener + tempro*enero)*deld*K.ergsev/2
        sumh2 += (tempr*(ener - eth) + tempro*(enero - eth))*deld*K.ergsev/2
        exptmp = expo(-max(fo_exptst_floor, (epii - eth)/bktm))
        bbnurj = min(epii, bb_energy_cap)^3*K.bb_coeff
        tempi1 = rnist*bbnurj*exptmp*sgtmp/epii
        tempi2 = rnist*bremtmp*exptmp*sgtmp/epii
        tempio = tempi
        tempi = tempi1 + tempi2
        atmp2o = atmp2
        atmp2 = tempi1*epii
        atmp22o = atmp22
        atmp22 = tempi1*(epii - eth)
        sumi += (tempi + tempio)*deld/2
        sumc += (atmp2 + atmp2o)*deld*K.ergsev/2
        sumc2 += (atmp22 + atmp22o)*deld*K.ergsev/2
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

# ---------------------------------------------------------------------------
# find53, milne and phint53hunt: the cross-section tables of the superlevels

const find53_floor = 1e-26           # table values below this are floored in the power-law extrapolation

"""
    interpolate_cross_section(ε, σ, e)

Cross section at the energy `e` (Ry above the threshold) of the table `σ(ε)` (XSTAR's `find53`):
linear between the points, a power law in the last interval, and zero outside the table.
"""
function interpolate_cross_section(ε::AbstractVector, σ::AbstractVector, e)
    n = length(ε)
    (n > 0 && e >= 0 && e <= ε[n]) || return 0.0
    j = clamp(searchsortedlast(ε, e), 1, n - 1)
    if j + 1 == n
        slope = log(max(σ[n], find53_floor)/max(σ[j], find53_floor))/
            log(max(ε[n], find53_floor)/max(ε[j], find53_floor))
        s = σ[j]*(e/ε[j])^slope
    else
        s = -σ[j]*(e - ε[j + 1])/(ε[j + 1] - ε[j]) + σ[j + 1]*(e - ε[j])/(ε[j + 1] - ε[j])
    end
    max(0.0, s)
end

const intin_exp_limit = 90.0
const intin_zero = 1e-3

# the two integrals ∫ x² e^{-x} and ∫ x³ e^{-x} of the Milne relation between s1 and s2 (XSTAR's `intin`)
function milne_terms(s1, s2, s0, T)
    K = constants()
    del = K.inv_kB/T
    s1, s2, s0 = s1*del, s2*del, s0*del
    ri2 = 0.0
    if s1 - s0 < intin_exp_limit
        ri2 = exp(s0 - s1)*((s1*s1 + 2s1 + 2) - exp(s1 - s2)*(s2*s2 + 2s2 + 2))/del/sqrt(del)
        s0 < intin_zero && s2 < intin_zero && s1 < intin_zero && (ri2 = 0.0)
    end
    rr = exp(s0 - s1)*(s1^3 - exp(s1 - s2)*s2^3)
    ri3 = (rr/del/sqrt(del) + 3ri2)/del
    (ri2, ri3)
end

const milne_tolerance = Float64(0.01f0)
const milne_flat = 1e-24

"""
    milne_recombination(T, ε, σ, E_th)

Recombination coefficient (cm³ s⁻¹) from the Milne relation for the cross-section table `σ(ε)`
(`ε` in Ry above the threshold `E_th`, in Ry) at the temperature `T` in K (XSTAR's `milne`). The
sum over the table stops once it changes by less than 1%.
"""
function milne_recombination(T, ε, σ, E_th)
    K = constants()
    nt = length(ε)
    st = (ε[1] + E_th)*K.Ry_erg
    sumo, s = 1.0, 0.0
    i = 1
    while abs(s - sumo) > milne_tolerance*s && i < nt
        i += 1
        s1 = (ε[i - 1] + E_th)*K.Ry_erg
        s2 = (ε[i] + E_th)*K.Ry_erg
        s2 < s1 && return 0.0                         # (XSTAR returns without setting the coefficient)
        v1, v2 = σ[i - 1], σ[i]
        if v1 != 0 || v2 != 0
            rb = (v2 - v1)/(s2 - s1 + milne_flat)
            ra = v2 - rb*s2
            ri2, ri3 = milne_terms(s1, s2, st, T)
            sumo = s
            s += ra*ri2 + rb*ri3
        end
    end
    s*K.milne_coeff
end

const hunt_tolerance = 0.01           # relative change of the integrals that ends the passes
const hunt_ln2 = Float64(0.69315f0)   # ln 2 as written in phint53hunt
const hunt_min_sum = 1e-24

"""
    photoionization_integrals_hunt(rad, E_th, ε, σ, T, swrat, xnx)

Photoionization and recombination integrals of a cross-section table `σ(ε)` (cm², `ε` in Ry above the
threshold `E_th` in eV) over the spectrum of `rad` (XSTAR's `phint53hunt`): the grid is sampled every
`2ᵏ` bins and refined by halving the step until the rates change by less than 1%. Returns `pirt`, `rrrt`,
`piht`, `rrcl`, `piht2`, `rrcl2` as `photoionization_integrals_fo`. Bins already computed in an earlier
pass are reused, except for the recombination weight measured from the threshold (`rrcl2`), which keeps the
value of the last bin computed in the pass, as in the original.
"""
function photoionization_integrals_hunt(rad::Radiation, E_th, ε, σ, T, swrat, xnx)
    K = constants()
    epi, bremsa = rad.E, rad.F
    ncn2 = length(epi)
    ntmp = length(ε)
    zero_sums = (; pirt=0.0, rrrt=0.0, piht=0.0, rrcl=0.0, piht2=0.0, rrcl2=0.0)
    ntmp <= 0 && return zero_sums
    eth = E_th
    nb1 = nbin(rad, eth) + 1
    numcon3 = ncn2 - max(grid_min_bins, ncn2 ÷ grid_guard_fraction)
    nb1 >= numcon3 && return zero_sums

    bktm = K.kB_cgs*T*T_unit/K.ergsev
    rnist = K.fo_saha*swrat/T/sqrt(T)

    # the span of the table, in a power of two bins
    nphint = nbin(rad, ε[ntmp]*K.Ry_eV_single + eth)
    ndelt = max(nphint - nb1, 1)
    itmp = trunc(Int, log(ndelt)/hunt_ln2 + 0.5)
    while true
        ndelt = 2^itmp
        nphint = nb1 + ndelt
        etst = nphint <= numcon3 ? (epi[nphint] - eth)/K.Ry_eV_single : 0.0
        if nphint > numcon3 || etst > ε[ntmp]
            itmp -= 1
            itmp > 1 && continue
        end
        break
    end
    nphint = min(nphint, ncn2)

    done = falses(ncn2)
    W = promote_type(eltype(σ), eltype(ε), typeof(T), typeof(swrat), typeof(xnx), eltype(bremsa))
    ansar1 = zeros(W, ncn2)
    ansar2 = zeros(W, ncn2)
    nskp = ndelt
    sumr = sumh = sumh2 = sumc = sumc2 = sumi = 0.0
    tst1 = tst2 = tst3 = tst4 = 0.0
    while (tst3 > hunt_tolerance || tst1 > hunt_tolerance || tst2 > hunt_tolerance ||
            tst4 > hunt_tolerance || sumi <= hunt_min_sum) && nskp > 1
        nskp = max(1, nskp ÷ 2)
        sumro, sumho, sumio, sumco = sumr, sumh, sumi, sumc
        sumr = sumh = sumh2 = sumc = sumc2 = sumi = 0.0
        tempr = tempi = atmp2 = atmp22 = 0.0
        ener = epi[nb1]
        kl = max(nb1 - 1, 1)
        while kl <= nphint
            enero = ener
            epii = epi[kl]
            ener = epii
            bremtmp = bremsa[kl]/K.eightpi
            tempio = tempi
            atmp2o = atmp2
            atmp22o = atmp22
            sgtmp = 0.0
            if ener >= eth
                if !done[kl]
                    sgtmp = interpolate_cross_section(ε, σ, (ener - eth)/K.Ry_eV_single)
                    ansar1[kl] = sgtmp
                    exptmp = expo(-(epii - eth)/bktm)
                    bbnurj = min(bb_energy_cap, epii)^3
                    tempi1 = rnist*bbnurj*sgtmp*exptmp*K.bb_coeff/epii
                    tempi2 = rnist*bremtmp*sgtmp*exptmp/epii
                    tempi = tempi1 + tempi2
                    atmp2 = tempi*epii
                    atmp22 = tempi*(epii - eth)
                    ansar2[kl] = atmp2
                else
                    sgtmp = ansar1[kl]
                    atmp2 = ansar2[kl]
                    tempi = atmp2/epii
                end
            end
            tempro = tempr
            tempr = K.eightpi*sgtmp*bremtmp/epii
            deld = ener - enero
            sumr += (tempr + tempro)*deld/2
            sumh += (tempr*ener + tempro*enero)*deld/2
            sumh2 += (tempr*(ener - eth) + tempro*(enero - eth))*deld/2
            sumi += (tempi + tempio)*deld/2
            sumc += (atmp2 + atmp2o)*deld/2
            sumc2 += (atmp22 + atmp22o)*deld/2
            done[kl] = true
            kl += nskp
        end
        tst3 = abs((sumio - sumi)/(sumio + sumi + hunt_min_sum))
        tst1 = abs((sumro - sumr)/(sumro + sumr + hunt_min_sum))
        tst2 = abs((sumho - sumh)/(sumho + sumh + hunt_min_sum))
        tst4 = abs((sumco - sumc)/(sumco + sumc + hunt_min_sum))
    end
    (; pirt=sumr, rrrt=xnx*sumi, piht=sumh*K.ergsev, rrcl=xnx*sumc*K.ergsev,
        piht2=sumh2*K.ergsev, rrcl2=xnx*sumc2*K.ergsev)
end

# ---------------------------------------------------------------------------
# superlevels (ucalc types 70 and 99)

const super_log_density_max = 8.0     # type 70 limits its second tabulated log₁₀ density to this
const super_density_cap = 1e8         # type 70 limits the density of a neutral ion to this
const super_xs_cap = 1e6              # type 70 caps the scaled cross section (Mb)
const super_xs_trim = Float64(1f-6)   # ... and drops the tail below this fraction of the first point
const rrrt_floor = 1e-48              # no rates if the recombination integral is below this

"""
    photoionize_superlevel(coef, cell, tabulate; radiation=NO_RADIATION, correct_energy=false, index=false)

Photoionization and recombination of a superlevel (XSTAR ucalc types 70 and 99). The cross section
`σ(ε)` is the table of the record scaled so that the Milne relation reproduces the tabulated
recombination coefficient at the density and temperature of `cell`; `tabulate(T, ntot, E_th)` returns
that coefficient `rec` (cm³ s⁻¹), `ε` (Ry above the threshold) and `σ` (Mb), for the temperature `T` in K, the
density `ntot`, and the threshold `E_th` in Ry. The cross section is then integrated over `radiation` with
`photoionization_integrals_hunt`, and the photoionization rate and the heating are scaled by the
ratio of the tabulated and integrated recombination rates. The recombination cooling is not scaled
(as in ucalc). With `correct_energy` the energies measured from the levels are scaled
too (type 99).

Returns `init` (the level, limited to the last real level `nlev` − 1), `final` (`nlev` + parent
level − 1), `frate` and `irate` (s⁻¹), `fenergy`, `ienergy`, `fenergy2`, `ienergy2` (erg s⁻¹) and
`opacity` (zero, the arrays are not filled). When the recombination integral is negligible ucalc returns
no rates but leaves the integrals `piht2` and `rrcl2` in its energy outputs; that is kept.
"""
function photoionize_superlevel(coef, cell::Cell, tabulate; radiation=NO_RADIATION,
    correct_energy=false, index=false)
    K = constants()
    levels = levels_of(coef)
    nlev = nlevels(levels, coef.ion)

    none = (; init=0, final=0, frate=0., irate=0., fenergy=0., ienergy=0.,
        fenergy2=0., ienergy2=0., opacity=0.)
    idest1 = min(Int(coef.level), nlev - 1)
    idest2 = max(nlev + Int(coef.parent.level) - 1, nlev)
    index && return (; none..., init=idest1, final=idest2)
    lo = get(levels, (coef.ion, idest1), nothing)
    cont = get(levels, (coef.ion, nlev), nothing)
    (lo === nothing || cont === nothing) && return none
    ggup = Float64(cont.g)
    e2 = Float64(cont.E)
    ett = abs(Float64(lo.E) - e2)
    if idest2 > nlev
        par = get(levels, (coef.parent.ion, Int(coef.parent.level)), nothing)
        par === nothing && return none
        ggup = Float64(par.g)
        e2 += Float64(par.E)                  # (ucalc reads the energy of this level from stale memory)
        ett = abs(Float64(lo.E) + Float64(par.E))
    end
    ggup <= min_g && return none
    swrat = Float64(lo.g)/ggup

    T = cell.T
    xnx = cell.nₑ
    rec, ε, xs = tabulate(T*T_unit, cell.ntot, ett/K.Ry_eV_coarse)   # (ucalc divides by 13.6)
    σ = max.(xs*Mb, 0.0)
    r = photoionization_integrals_hunt(radiation, ett, ε, σ, T, swrat, xnx)
    ans = if r.rrrt <= rrrt_floor
        (0.0, 0.0, 0.0, 0.0, r.piht2, r.rrcl2)
    else
        scale = rec*xnx/r.rrrt
        a1, a2 = r.pirt*scale, rec*xnx
        a3, a4, a5, a6 = -r.rrcl, -r.piht*scale, -r.rrcl2, -r.piht2*scale
        if correct_energy
            dE = (e2 - Float64(lo.E))*K.ergsev
            a6 *= (abs(a4) - dE*a1)/max(heat_floor, abs(a4) - ett*K.ergsev*a1)
            a5 *= (abs(a3) - dE*a2)/max(heat_floor, abs(a3) - ett*K.ergsev*a2)
        end
        (a1, a2, a3, a4, a5, a6)
    end
    (; init=idest1, final=idest2, frate=ans[1], irate=ans[2], fenergy=-ans[4],
        ienergy=-ans[3], fenergy2=-ans[6], ienergy2=-ans[5], opacity=0.0)
end

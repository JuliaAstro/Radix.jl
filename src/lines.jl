# The opacity of a line put into the continuum bins: XSTAR's linopac, with the Voigt function voigte it uses.

const voigt_a_zero = 0.0
const voigt_a_small = Float64(0.2f0)
const voigt_a_large = Float64(1.4f0)
const voigt_u_large = Float64(3.2f0)
const voigt_v_far = 5.0
const voigt_v_exp = 100.0               # exp(-v²) is taken to be 0 beyond this v²
const voigt_tail_shift = Float64(1.5f0)
const voigt_v_mid = Float64(2.4f0)
const voigt_v_low = Float64(1.3f0)
const voigt_sqrt_pi = Float64(1.772453851f0)
const voigt_sqrt_2 = Float64(1.414213562f0)
# the coefficients are single-precision literals in XSTAR's source
const voigt_ak = Float64.((-1.12470432f0, -0.15516677f0, 3.28867591f0, -2.34357915f0, 0.42139162f0, -4.48480194f0, 9.39456063f0, -6.61487486f0, 1.98919585f0, -0.22041650f0, 0.554153432f0, 0.278711796f0, -0.188325687f0, 0.042991293f0, -0.003278278f0, 0.979895023f0, -0.962846325f0, 0.532770573f0, -0.122727278f0))

"""
    voigt(v, a)

The Voigt function `H(a, v)` of the line profile with the Doppler width as the unit of `v` and the damping parameter `a`, XSTAR's `voigte`
(the approximation of Hui, Armstrong and Wray in the pieces that depend on `a` and `v`; `exp(-v²)` for `a = 0`).
"""
function voigt(vs, a)
    v = abs(vs)
    u = a + v
    v2 = v*v
    a == voigt_a_zero && return v2 < voigt_v_exp ? exp(-v2) : zero(v2)
    if a <= voigt_a_small
        v >= voigt_v_far && return a*(15 + 6v2 + 4v2*v2)/(4v2*v2*v2*voigt_sqrt_pi)
        return voigt_poly(v, v2, a, v2 < voigt_v_exp ? exp(-v2) : zero(v2), 1)
    elseif a > voigt_a_large || u > voigt_u_large
        a2 = a*a
        u = voigt_sqrt_2*(a2 + v2)
        u2 = 1/(u*u)
        return voigt_sqrt_2/voigt_sqrt_pi*a/u*(1 + u2*(3v2 - a2) + u2*u2*(15v2*v2 - 30v2*a2 + 3a2*a2))
    else
        return voigt_poly(v, v2, a, v2 < voigt_v_exp ? exp(-v2) : zero(v2), 2)
    end
end

# the piece of voigte for a ≤ 0.2 (k = 1) or 0.2 < a ≤ 1.4 with a + v ≤ 3.2 (k = 2), with the polynomial of the range of v
function voigt_poly(v, v2, a, ex, k)
    ak = voigt_ak
    quo = one(v)
    m = v < voigt_v_mid ? (v < voigt_v_low ? 1 : 6) : 11
    v >= voigt_v_mid && (quo = 1/(v2 - voigt_tail_shift))
    h1 = quo*(ak[m] + v*(ak[m + 1] + v*(ak[m + 2] + v*(ak[m + 3] + v*ak[m + 4]))))
    k == 1 && return h1*a + ex*(1 + a*a*(1 - 2v2))
    pqs = 2/voigt_sqrt_pi
    h1p = h1 + pqs*ex
    h2p = pqs*h1p - 2v2*ex
    h3p = (pqs*(1 - ex*(1 - 2v2)) - 2v2*h1p)/3 + pqs*h2p
    h4p = (2v2*v2*ex - pqs*h1p)/3 + pqs*h3p
    psi = ak[16] + a*(ak[17] + a*(ak[18] + a*ak[19]))
    psi*(ex + a*(h1p + a*(h2p + a*(h3p + a*h4p))))
end

# ---------------------------------------------------------------------------
# linopac

const line_grid_size = 20000            # nbtpp: points of the temporary grid on which the profile is calculated
const line_bins_max = 999999            # ncn: the size of XSTAR's arrays of the energy grid
const line_wavelength_max = 1e8         # Å; lines outside this range and 1 Å are not put in the bins
const line_wavelength_min = 1.0
const line_hc_floor = 1e-49
const line_damping_floor = Float64(1f-24)         # aasmall: added to the width
const line_profile_norm = Float64(1.772f0)        # √π as linopac writes it
const line_damping_voigt = Float64(1f-6)          # aasmall above which the profile at the bin of the line centre is a Voigt function
const line_damping_voigt_wing = Float64(1f-9)     # the same on the temporary grid
const line_profile_cut = Float64(1f-6)            # dpcrit: the profile is not followed below this value
const line_width_max = 50.0                       # in Doppler widths
const line_width_damping = Float64(200f0)         # or in this many damping parameters
const line_sum_floor = 1e-34
const line_single_bin = 2                         # lfast above this puts the line in one bin

struct LineBuffers
    energy::Vector{Float64}
    opacity::Vector{Float64}
end

line_buffers() = get!(() -> LineBuffers(zeros(line_grid_size), zeros(line_grid_size)), task_local_storage(), :RadixLineBuffers)::LineBuffers

"""
    add_line!(opacity, radiation, optpp, rcem1, rcem2, elin, vturb, T, mass, delea; lfast=2)

Puts the line of wavelength `elin` (Å) and centre opacity `optpp` (cm⁻¹, `opacity` of the rate of the line times the density of its lower
level) into the bins of `radiation.E` in the `opacity` (`Opacity`: `total`), as XSTAR's `linopac` does. `vturb` is the turbulent speed (km/s), `T`
the temperature (10⁴ K), `mass` the mass of the ion (amu) and `delea` the natural width (eV). The profile is the Voigt function of the
thermal and turbulent widths and the damping parameter `delea/(width 4π)`, calculated on a grid `ncut` times finer than the bin of the line
centre and averaged in each bin (`lfast ≤ 2`); with `lfast > 2` all the opacity goes in the bin of the centre, in `total` as `optpp × width` over
the width of the bin, and the emissivities `rcem1` and `rcem2` (erg s⁻¹ cm⁻³) in `emissivity` (they are not used for the full profile).
"""
function add_line!(opacity::Opacity, radiation::Radiation, optpp, rcem1, rcem2, elin, vturb, T, mass, delea; lfast=2)
    K = constants()
    E = radiation.E
    ncn2 = length(E)
    (elin > line_wavelength_max || elin < line_wavelength_min) && return opacity
    elin = abs(elin)
    vth = K.thermal_speed/cm_per_km*sqrt(T/mass)
    e0 = K.hc_eVÅ_single/max(elin, line_hc_floor)
    e0 <= E[1] && return opacity
    deleturb = e0*(vturb/(K.light_speed/cm_per_km))
    deleth = e0*(vth/(K.light_speed/cm_per_km))
    dele = sqrt(deleth*deleth + deleturb*deleturb)
    aasmall = delea/(line_damping_floor + dele)/K.fourpi

    ml1 = max(min(line_bins_max - 1, nbin(radiation, e0)), 2)
    prftmp = 2/(E[ml1 + 1] - E[ml1 - 1])
    opsv4 = optpp*dele
    if lfast > line_single_bin
        opacity.total[ml1] += opsv4*prftmp
        opacity.emissivity[1, ml1] += rcem1*prftmp/K.ergsev/K.fourpi
        opacity.emissivity[2, ml1] += rcem2*prftmp/K.ergsev/K.fourpi
        return opacity
    end

    buffers = line_buffers()
    etpp, optpp2 = buffers.energy, buffers.opacity
    e00 = E[ml1]
    deleepi = E[ml1 + 1] - E[ml1]
    ncut = min(max(trunc(Int, deleepi/dele), 1), line_grid_size ÷ 10)
    deleused = deleepi/ncut
    ml2 = line_grid_size ÷ 2
    mlmin, mlmax = line_grid_size, 1
    ml1min, ml1max = line_bins_max + 1, 0
    done = (false, false)
    delet = (e00 - e0)/dele
    profile = (aasmall > line_damping_voigt ? voigt(abs(delet), aasmall) : exp(-delet*delet))/line_profile_norm
    etpp[ml2] = e00
    optpp2[ml2] = optpp*profile
    tst = one(profile)
    mlc = 0
    ldir = 1
    while !(done[1] && done[2]) && mlc < line_grid_size ÷ 2          # (the test of the end of the loop is never met: ml1min stays above ml1)
        mlc += 1
        for ij in 1:2
            ldir = -ldir
            done[ij] && continue
            mlm = ml2 + ldir*mlc
            etptst = e00 + ldir*mlc*deleused
            if mlm <= line_grid_size && mlm >= 1 && etptst > 0 && etptst < E[ncn2]
                mlmin = min(mlm, mlmin)
                mlmax = max(mlm, mlmax)
                etpp[mlm] = e00 + ldir*mlc*deleused
                delet = (etpp[mlm] - e0)/dele
                profile = (aasmall > line_damping_voigt_wing ? voigt(abs(delet), aasmall) : exp(-delet*delet))/line_profile_norm
                optpp2[mlm] = optpp*profile
                tst = profile
            end
            if (tst < line_profile_cut || mlm <= 1 || mlm >= line_grid_size || etptst <= 0 || etptst >= E[ncn2] ||
                    mlc > line_grid_size || abs(delet) > max(line_width_max, line_width_damping*aasmall)) &&
                    ml1min < ml1 - 2 && ml1max > ml1 + 2 && ml1min >= 1 && ml1max <= line_bins_max
                done = ij == 1 ? (true, done[2]) : (done[1], true)
            end
        end
    end

    # average the profile in the bins of the energy grid
    opsum = tmpop = sume = 0.0
    ml1m = nbin(radiation, etpp[mlmin])
    mlmin = max(mlmin, 2)
    mlmax = min(mlmax, line_grid_size)
    for mlm in mlmin + 1:mlmax
        tmpopo = tmpop
        tmpop = optpp2[mlm]
        tmpe = abs(etpp[mlm] - etpp[mlm - 1])
        sume += tmpe
        opsum += (tmpop + tmpopo)*tmpe/2
        if etpp[mlm] > E[ml1m]
            if sume > line_sum_floor
                average = opsum/sume
                while etpp[mlm] > E[ml1m] && ml1m < ncn2
                    opacity.total[ml1m] += average
                    ml1m += 1
                end
            end
            opsum = sume = 0.0
        end
    end
    opacity
end

# ---------------------------------------------------------------------------
# which lines go in the continuum bins, as `ucalc` and `calc_emis_ion` decide

const line_data_type = 4                          # the records of this type have an emissivity and are ranked by it (`calc_emisab_ion`)
const line_opacity_min = Float64(1f-34)           # ucalc (type 50): the lines weaker than this are not put in the bins
const line_abundance_min = Float64(1f-34)         # calc_emisab_ion: lines whose levels are both emptier than this have no emissivity
const line_bin_floor = 1e-34                      # rlbin: added to the wavelength of the line for its bin
const line_rank_floor = Float64(1f-37)            # rlbin: lines with a centre opacity and an emissivity below this are not ranked
const line_width_per_A = Float64(4.136f-15)       # ucalc (type 50): the natural width (eV) is the Einstein A (s⁻¹) times this (h in eV s)
const default_turbulence = 1.0                    # km/s, as in `rate` of the lines

"""
    rank!(list, entry, nrank)

Puts `entry = (strength, k, j)` in the `list` of the strongest entries of an energy bin, as XSTAR's `rlbin` does: the list is in order of
decreasing strength and holds `nrank` entries at most. An entry whose place would be the last one (`nrank`) is not inserted, as in
`rlbin`, which returns before it stores it; the entries below it move down and the last one drops out.
"""
function rank!(list::AbstractVector, entry, nrank)
    mm = 0
    while true
        mm += 1
        mm > length(list) && break
        entry[1] > list[mm][1] && break
        mm >= nrank && break
    end
    mm >= nrank && return list
    insert!(list, mm, entry)
    length(list) > nrank && pop!(list)
    list
end

# the natural width (eV) with which the record puts its line in the bins, or nothing for the records that do not (only type 50 does)
line_width(::AbstractRate, widths=nothing) = nothing
function line_width(coef::AtomicLine2, widths=nothing)
    width = Float64(coef.A)*line_width_per_A
    widths === nothing && return width
    a, b = get(coef.levels, (coef.ion, coef.transition.lower), nothing), get(coef.levels, (coef.ion, coef.transition.upper), nothing)
    (a === nothing || b === nothing) && return width
    upper = a.E < b.E ? b.level : a.level                              # (the level `idest1` of ucalc, the higher by energy)
    auger = get(widths, (Int(coef.ion), Int(upper)), nothing)         # deleafnd: the width of a K-vacancy level of the database
    auger === nothing ? width : auger*line_width_per_A
end

"""
    auger_widths(mixture)

The widths `A_auto(k, parent)` (s⁻¹, the third real of the record) of the K-vacancy levels of the records of rate type 41 (`IronKAuger`) of the elements of `mixture`, as a dictionary `(ion, level) => width` which `deleafnd` consults
for the natural width of a line whose upper level is one of them (the first record of the level, as it does).
"""
function auger_widths(mixture)
    Dict{Tuple{Int, Int}, Float64}(key => first(rates) for (key, rates) in auger_rates(mixture))
end

"""
    auger_rates(mixture)

The rates `(A_auto(k, parent), A_rad(k))` (s⁻¹, the third and fourth reals) of the first record of rate type 41 of each K-vacancy level `(ion, level)`: `auger_widths` has the first, `binemis` takes both.
"""
function auger_rates(mixture)
    rates = Dict{Tuple{Int, Int}, Tuple{Float64, Float64}}()
    for element in mixture.elements, coef in element.rates
        coef isa IronKAuger || continue
        get!(rates, (Int(coef.ion), Int(coef.level)), (Float64(coef.A_widths[2]), Float64(coef.A_widths[3])))
    end
    rates
end
has_line(::AbstractRate) = false
has_line(::Union{AtomicLine2, RadiativeAPED}) = true

"""
    add_line!(opacity, radiation, coef, centre, emissivity, T; vturb=default_turbulence, lfast=2)

Puts the line of the record `coef` in the bins with `add_line!`, as `ucalc` does for a line of type 50 (the other records do nothing): `centre` is the opacity at the
line centre (`oplin`, cm⁻¹) and `emissivity` the pair of emissivities `(rcem1, rcem2)`. The natural width is the Einstein A times 4.136×10⁻¹⁵ eV s, or, for a line whose upper level is a K-vacancy level with a record of rate type 41, the width of that record (`widths` of `auger_widths`: `deleafnd`).
"""
function add_line!(opacity::Opacity, radiation::Radiation, coef::AtomicLine2, centre, emissivity, T; vturb=default_turbulence, lfast=2, widths=nothing)
    centre > line_opacity_min || return opacity
    mass = atomic_mass(coef.levels, coef.ion)
    add_line!(opacity, radiation, centre, emissivity[1], emissivity[2], abs(Float64(coef.λ)), vturb, T, mass, line_width(coef, widths); lfast)
end

const line_fe_opacity_min = 1e-48                 # ucalc (type 82): the lines weaker than this are not put in the bins

# the Fe UTA lines (type 82): `ucalc` calls linopac with no emission and the width `rdat1(np1r-1+6)*4.14e-15` eV, the sixth real of a record that has five (the first of the next record, a wavelength: about 10⁻¹³ eV, much
# less than the Doppler width, so that the profile is a Gaussian, which is what a width of 0 gives)
function add_line!(opacity::Opacity, radiation::Radiation, coef::RadiativeFeDecay, centre, emissivity, T; vturb=default_turbulence, lfast=2, widths=nothing)
    centre > line_fe_opacity_min || return opacity
    mass = atomic_mass(coef.levels, coef.ion)
    add_line!(opacity, radiation, centre, 0.0, 0.0, abs(Float64(coef.λ)), vturb, T, mass, 0.0; lfast)
end
add_line!(opacity::Opacity, radiation::Radiation, ::AbstractRate, centre, emissivity, T; kw...) = opacity

# ---------------------------------------------------------------------------
# two-photon decays

const two_photon_data_type = 9                    # the lines of this type are two-photon decays: they emit a continuum, not a line
const two_photon_size_min = Float64(0.01f0)       # calc_emis_ion: the records of type 9 whose first real (the wavelength of a line, the rate of a decay) is below this are not called
two_photon_size(coef::AtomicLine2) = Float64(coef.λ)
two_photon_size(coef::TwoPhotonDecay) = Float64(coef.A)
const two_photon_sum_floor = 1e-24                # ucalc: added to the integral of the shape in its normalization

"""
    add_two_photon!(opacity, radiation, coef, abund2, ptmp)

Puts the continuum of the two-photon decay `coef` (a `TwoPhotonDecay`, or an `AtomicLine2` of rate type 9) in the recombination emissivity of `opacity`, as `ucalc` does (`ind = 50`, `nrdesc = 9`):
the energy `A hν` of the decay of the upper level, with the density `abund2`, is spread over the bins below the energy `hν` of the transition with the shape `E² (hν - E)`,
`ptmp[1]` of it in the first direction and `ptmp[2]` in the second (per steradian: the emissivity is `abund2 A hν shape ptmp / (4π ∫ shape dE)`; the integral starts from 0 at the first bin as `ucalc` does it). The records with a level out of the table
or a wavelength (or rate) below 0.01 add nothing.
"""
function add_two_photon!(opacity::Opacity, radiation::Radiation, coef::Union{AtomicLine2, TwoPhotonDecay}, abund2, ptmp)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    i1, i2 = coef.transition.lower, coef.transition.upper
    (i1 <= 0 || i1 >= nlev || i2 <= 0 || i2 >= nlev || two_photon_size(coef) <= two_photon_size_min) && return opacity
    a, b = get(levels, (coef.ion, i1), nothing), get(levels, (coef.ion, i2), nothing)
    (a === nothing || b === nothing) && return opacity
    K = constants()
    E = radiation.E
    emax = abs(Float64(a.E) - Float64(b.E))
    nbmx = nbin(radiation, emax)
    shape(ll) = E[ll]^2*max(0.0, E[nbmx] - E[ll])
    total = 0.0
    previous = 0.0
    for ll in 2:nbmx
        current = shape(ll)
        total += (current + previous)*(E[ll] - E[ll - 1])/2
        previous = current
    end
    scale = Float64(coef.A)*emax/(two_photon_sum_floor + total)
    for ll in 2:nbmx
        emission = abund2*shape(ll)*scale/K.fourpi
        opacity.emissivity[1, ll] += emission*ptmp[1]
        opacity.emissivity[2, ll] += emission*ptmp[2]
    end
    opacity
end

# ---------------------------------------------------------------------------
# the emission of the lines binned in the spectrum: binemis

const broaden_velocity = Float64(1.2f1)         # binemis: the thermal speed of an atom of 1 amu at 10⁴ K (km/s; linopac has 12.9)
const broaden_hc = Float64(12398.42f0)          # binemis: the energy of the wavelength of a line
const broaden_light = 3e5                       # km/s
const broaden_width = Float64(4.14f-15)         # the widths in eV of the rates in s⁻¹
const broaden_floor = Float64(1f-36)            # added to the Doppler width in the damping parameter
const broaden_unit = Float64(1.602197f-12)      # the profile is per erg
const broaden_voigt_min = Float64(1f-9)         # damping parameter above which the profile is a Voigt function
const broaden_sum_floor = 1e-24                 # the interval of the average of a bin is at least this (eV)
const broaden_buffer = 999999                   # ncn: the points of the profile are in an array of this size
const broaden_bin_floor = Float64(1f-34)        # added to the wavelength for the energy of the line

"""
    broaden(E, lines, T, vturb)

The emission of lines on the energy grid `E` (eV), as XSTAR's `binemis` makes it: `lines` is a collection of `(wavelength, mass, inward, outward, delea, egam)` (Å, amu, the two luminosities in erg s⁻¹, the Auger width in s⁻¹ and the radiative width
`egam` that sets the damping), `T` the temperature (10⁴ K) and `vturb` the turbulent speed (km s⁻¹; the larger of it and the thermal speed is used). Each line has a Voigt profile, calculated on a grid `ncut` times finer than the bin of the line
and averaged in the bins of `E`, per erg. Returns `(inward, outward)`, in erg s⁻¹ erg⁻¹.
"""
function broaden(E::AbstractVector, lines, T, vturb)
    K = constants()
    n = length(E)
    rad = Radiation(E, zeros(n))
    inward, outward = zeros(n), zeros(n)
    buffer = broaden_buffer
    z1, z2 = zeros(buffer), zeros(buffer)
    s1, s2 = zeros(n), zeros(n)               # (zrtmps: set by the bins of a line, kept for the next one as `binemis` does)
    for (elin, mass, lum1, lum2, delea, egam) in lines
        eline = K.hc_eVÅ_single/(broaden_bin_floor + elin)
        nb1 = nbin(rad, eline)
        nb1 > 2 || continue
        vth = broaden_velocity*sqrt(T/mass)
        vt = max(vturb, vth)
        e0 = broaden_hc/max(elin, line_hc_floor)
        dele = hypot(e0*(vt/broaden_light), e0*(vth/broaden_light))
        aasmall = (delea + egam*broaden_width)/(broaden_floor + dele)/K.fourpi
        profile(x) = (aasmall > broaden_voigt_min ? voigt(abs(x), aasmall) : exp(-x*x))/line_profile_norm/dele/broaden_unit
        e00 = E[nb1]
        deleepi = E[nb1 + 1] - E[nb1]
        ncut = min(max(trunc(Int, deleepi/dele), 1), buffer ÷ 10)
        deleused = deleepi/ncut
        ml2 = buffer ÷ 2
        p = profile((e00 - e0)/dele)
        z1[ml2], z2[ml2] = lum1*p, lum2*p
        # the points on the two sides of the centre, `binemis` goes on to 249999 steps and does nothing outside the grid
        mlmin, mlmax = buffer, 1
        top = E[n]
        @inbounds for mlc in 1:buffer ÷ 2 - 1
            etptst = e00 - mlc*deleused
            etptst > 0 || break
            mlm = ml2 - mlc
            mlm > 1 || break
            mlm < mlmin && (mlmin = mlm)
            mlm > mlmax && (mlmax = mlm)
            p = profile((etptst - e0)/dele)
            z1[mlm], z2[mlm] = lum1*p, lum2*p
        end
        @inbounds for mlc in 1:buffer ÷ 2 - 1
            etptst = e00 + mlc*deleused
            etptst < top || break
            mlm = ml2 + mlc
            mlm < buffer || break
            mlm < mlmin && (mlmin = mlm)
            mlm > mlmax && (mlmax = mlm)
            p = profile((etptst - e0)/dele)
            z1[mlm], z2[mlm] = lum1*p, lum2*p
        end
        mlmin <= mlmax || continue
        point(m) = e00 + (m - ml2)*deleused           # (the energy of a point of the profile)
        ml1min, ml1max = nbin(rad, point(mlmin)), nbin(rad, point(mlmax))
        ml1m = ml1min
        mlmin, mlmax = max(mlmin, 2), min(mlmax, buffer)
        sume = sum1 = sum2 = 0.0
        @inbounds for mlm in mlmin + 1:mlmax
            xm = point(mlm)
            tmpe = abs(xm - point(mlm - 1))
            sume += tmpe
            sum1 += (z1[mlm] + z1[mlm - 1])*tmpe/2
            sum2 += (z2[mlm] + z2[mlm - 1])*tmpe/2
            if xm > E[ml1m]
                mlm == mlmax && (ml1m = max(1, ml1m - 1))
                if sume > broaden_sum_floor
                    a1, a2 = sum1/sume, sum2/sume
                    while xm > E[ml1m] && ml1m < n
                        s1[ml1m], s2[ml1m] = a1, a2
                        ml1m += 1
                    end
                end
                sum1 = sum2 = sume = 0.0
            end
        end
        z1[mlmin:mlmax] .= 0
        z2[mlmin:mlmax] .= 0
        for m in ml1min:ml1max
            inward[m] += s1[m]
            outward[m] += s2[m]
        end
        for m in ml1min:ml1max
            s1[m] = s2[m] = 0.0
        end
    end
    (; inward, outward)
end

"""
    broaden(model, run; vturb=run.vturb, min_fraction=1e-15, nrank=rank_per_bin)

The emission of the lines of a run (`slab_model`) binned on its energy grid (`broaden(E, lines, T, vturb)`), at the temperature of the last zone: the `nrank` lines with the largest luminosity in each bin
(`rank!`) among those with more than `min_fraction` of the luminosity of the source. Returns `(inward, outward)` in erg s⁻¹ erg⁻¹.
"""
function broaden(model, run; vturb=run.vturb, min_fraction=results_min_fraction, nrank=rank_per_bin)
    mixture = model.mixture
    E = run.E
    rad = Radiation(E, zeros(length(E)))
    threshold = min_fraction*band_luminosity(E, run.incident)
    rates = auger_rates(mixture)
    emina, emaxa = constants().hc_eVÅ_single/E[end], constants().hc_eVÅ_single/E[1]
    ranked = Dict{Int, Vector{Tuple{Float64, Int, Int}}}()
    for (k, element) in enumerate(mixture.elements), j in eachindex(element.rates)
        c = element.rates[j]
        (c isa AtomicLine2 && c.rtype == line_data_type && element.lo[j] != 0) || continue
        inward, outward = run.luminosities.inward[k][j]*luminosity_unit, run.luminosities.outward[k][j]*luminosity_unit
        (inward > threshold || outward > threshold) || continue
        elin = abs(Float64(c.λ))
        (emina <= elin <= emaxa) || continue
        emission = (inward + outward)/luminosity_unit
        emission < line_rank_floor && continue
        rank!(get!(ranked, nbin(rad, constants().hc_eVÅ_single/(line_bin_floor + elin)), Tuple{Float64, Int, Int}[]), (emission, k, j), nrank)
    end
    selected = NTuple{6, Float64}[]               # (a concrete type: the loop of `broaden` is type-stable and fast)
    for bin in sort!(collect(keys(ranked))), (_, k, j) in ranked[bin]
        c = mixture.elements[k].rates[j]
        auger = get(rates, (Int(c.ion), Int(c.transition.upper)), nothing)
        delea, egam = auger === nothing ? (0.0, Float64(c.A)) : (auger[1]*broaden_width, auger[2])
        push!(selected, (abs(Float64(c.λ)), Float64(atomic_mass(c.levels, c.ion)), run.luminosities.inward[k][j]*luminosity_unit,
            run.luminosities.outward[k][j]*luminosity_unit, delea, egam))
    end
    broaden(E, selected, run.zones[end].T, vturb)
end

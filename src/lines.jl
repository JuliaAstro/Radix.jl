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
    ex = v2 < voigt_v_exp ? exp(-v2) : zero(v2)
    if a <= voigt_a_small
        v >= voigt_v_far && return a*(15 + 6v2 + 4v2*v2)/(4v2*v2*v2*voigt_sqrt_pi)
        return voigt_poly(v, v2, a, ex, 1)
    elseif a > voigt_a_large || u > voigt_u_large
        a2 = a*a
        u = voigt_sqrt_2*(a2 + v2)
        u2 = 1/(u*u)
        return voigt_sqrt_2/voigt_sqrt_pi*a/u*(1 + u2*(3v2 - a2) + u2*u2*(15v2*v2 - 30v2*a2 + 3a2*a2))
    else
        return voigt_poly(v, v2, a, ex, 2)
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
line_width(::AbstractRate) = nothing
line_width(coef::AtomicLine2) = Float64(coef.A)*line_width_per_A
has_line(::AbstractRate) = false
has_line(::Union{AtomicLine2, RadiativeAPED}) = true

"""
    add_line!(opacity, radiation, coef, centre, emissivity, T; vturb=default_turbulence, lfast=2)

Puts the line of the record `coef` in the bins with `add_line!`, as `ucalc` does for a line of type 50 (the other records do nothing): `centre` is the opacity at the
line centre (`oplin`, cm⁻¹) and `emissivity` the pair of emissivities `(rcem1, rcem2)`. The natural width is the Einstein A times 4.136×10⁻¹⁵ eV s: XSTAR
takes it from the records of type 41 of the database where there are some, which are not read here.
"""
function add_line!(opacity::Opacity, radiation::Radiation, coef::AtomicLine2, centre, emissivity, T; vturb=default_turbulence, lfast=2)
    centre > line_opacity_min || return opacity
    mass = atomic_mass(coef.levels, coef.ion)
    add_line!(opacity, radiation, centre, emissivity[1], emissivity[2], abs(Float64(coef.λ)), vturb, T, mass, line_width(coef); lfast)
end
add_line!(opacity::Opacity, radiation::Radiation, ::AbstractRate, centre, emissivity, T; kw...) = opacity

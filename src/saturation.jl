# The absorption of the lines in the continuum bins: the lines are narrower than a bin at any number of bins (the thermal width is about 6×10⁻⁵ of the energy), so that the
# opacity averaged over a bin, which is what XSTAR's `linopac` puts in the arrays, makes the transmission `exp(-⟨τ⟩)` instead of `⟨exp(-τ)⟩` and overestimates the absorption of
# a saturated line (the equivalent width of the slab with its lines 1.5 to 3 times that of 30000 bins at 999). `LineAbsorbers` keeps the strong lines of the zones with their profiles
# (the integral of the opacity, the widths) and `resolved_transmission` gives the transmission of the bins with the equivalent width of each line in the cells of the bins.

const absorber_peak_min = 1e-6          # the lines with a centre depth below this are not followed
const absorber_thin = 0.05              # the lines with a centre depth below this put all their mean depth `S dele` in the cell of their centre
const absorber_saturated = 1.0          # the lines with a centre depth above this are treated with their equivalent width, the others (up to it) with their mean depth in each cell
const absorber_step_core = 0.25         # the profile is followed in steps of this many Doppler widths up to absorber_core
const absorber_core = 4.0
const absorber_step_wing = 1.05         # and then with steps that grow by this factor

"""
    LineAbsorbers()

The strong lines of the zones of a slab, summed over the zones (`add_zone!`): for each line the centre energy `e0`, the integral `S` over the zones of the centre opacity times the thickness (`optpp Δr`, so that the depth of the
line is `S H(a, v)/√π` with the Voigt function `H`), and the mean Doppler width and damping parameter, weighted by `S`.
"""
struct LineAbsorbers
    index::Dict{Tuple{Int, Int}, Int}
    e0::Vector{Float64}
    S::Vector{Float64}
    widths::Vector{Float64}          # Σ S dele
    damping::Vector{Float64}         # Σ S a
end
LineAbsorbers() = LineAbsorbers(Dict{Tuple{Int, Int}, Int}(), Float64[], Float64[], Float64[], Float64[])

# the energy, the Doppler width (eV) and the damping parameter with which `add_line!` puts the line of the record in the bins
function line_profile(coef::Union{AtomicLine2, RadiativeFeDecay}, T, vturb, widths)
    K = constants()
    elin = abs(Float64(coef.λ))
    (elin > line_wavelength_max || elin < line_wavelength_min) && return nothing
    mass = weight(coef.levels, coef.ion)
    delea = coef isa AtomicLine2 ? line_width(coef, widths) : 0.0
    e0 = K.hc_eVÅ_single/max(elin, line_hc_floor)
    c = K.light_speed/cm_per_km
    dele = sqrt((e0*K.thermal_speed/cm_per_km*sqrt(T/mass)/c)^2 + (e0*vturb/c)^2)
    (e0, dele, delea/(line_damping_floor + dele)/K.fourpi)
end
line_profile(::AbstractRate, T, vturb, widths) = nothing

"""
    add_zone!(absorbers, continuum, edges, T, Δr)

Adds the lines of a zone of thickness `Δr` (cm) and temperature `T` whose centre opacities are `edges` (`record_opacities`) to `absorbers`: the records of `continuum` (`Continuum`) that `add_line!` puts in
the bins, whose centre depth `optpp Δr/√π` is above `absorber_peak_min`.
"""
function add_zone!(absorbers::LineAbsorbers, continuum::Continuum, edges, T, Δr)
    for (k, j, _) in continuum.catalog.lines
        centre = edges[k][j]
        centre*Δr/line_profile_norm > absorber_peak_min || continue
        coef = continuum.mixture.elements[k].rates[j]
        profile = line_profile(coef, T, continuum.vturb, continuum.widths)
        profile === nothing && continue
        e0, dele, a = profile
        i = get!(absorbers.index, (k, j)) do
            push!(absorbers.e0, e0); push!(absorbers.S, 0.0); push!(absorbers.widths, 0.0); push!(absorbers.damping, 0.0)
            length(absorbers.e0)
        end
        S = centre*Δr
        absorbers.S[i] += S
        absorbers.widths[i] += S*dele
        absorbers.damping[i] += S*a
    end
    absorbers
end

"""
    resolved_transmission(E, absorbers, depth)

The transmission of the bins of the grid `E` of a slab: `depth` is the depth of the slab without the lines (`dpthc` minus the mean depth that the lines put in the bins, the opacity averaged in each bin) and the lines
are those of `absorbers` with their profiles inside the cells of the bins (from the middle of the interval to the bin before to the middle of that to the one after, the weights of the trapezoid rule). A line with a centre depth `S/√π`
below `absorber_thin` puts its mean depth `S dele` all in the cell of its centre; one above `absorber_saturated` removes the fraction `W/Δ` of the cell, `W = ∫ (1 - e^{-τ}) dE` of its profile in the cell of width `Δ`, and the others add their
mean depth in each cell: the bin is `exp(-(depth + Σ ⟨τ⟩)) (1 - Σ W/Δ)`. This is the average of `exp(-τ)` over the bin, and not `exp(-⟨τ⟩)`, which overestimates the absorption of a saturated line;
overlapping saturated lines in a cell are not combined (the transmission of the cell is at least 0).
"""
function resolved_transmission(E::AbstractVector, absorbers::LineAbsorbers, depth::AbstractVector)
    n = length(E)
    mean_depth = zeros(n)
    removed = zeros(n)
    cell_lo(k) = k == 1 ? E[1] : (E[k - 1] + E[k])/2
    cell_hi(k) = k == n ? E[n] : (E[k] + E[k + 1])/2
    cell_of(x) = (k = clamp(searchsortedlast(E, x), 1, n - 1); x >= (E[k] + E[k + 1])/2 ? k + 1 : k)
    for i in eachindex(absorbers.e0)
        S = absorbers.S[i]
        e0, dele, a = absorbers.e0[i], absorbers.widths[i]/S, absorbers.damping[i]/S
        (E[1] < e0 < E[n]) || continue
        peak = S*voigt(0.0, a)/line_profile_norm
        if peak <= absorber_thin
            k = cell_of(e0)
            mean_depth[k] += S*dele/(cell_hi(k) - cell_lo(k))
            continue
        end
        saturated = peak > absorber_saturated
        xmax = max(line_width_max, line_width_damping*a)
        τ(x) = S*voigt(x, a)/line_profile_norm
        x, τx, h = 0.0, τ(0.0), absorber_step_core
        gx = -expm1(-τx)
        while x < xmax
            xn = x + h
            τn = τ(xn)
            gn = -expm1(-τn)
            for side in (-1, 1)
                k = cell_of(e0 + side*dele*(x + xn)/2)
                width = cell_hi(k) - cell_lo(k)
                if saturated
                    removed[k] += dele*(gx + gn)/2*(xn - x)/width
                else
                    mean_depth[k] += dele*(τx + τn)/2*(xn - x)/width
                end
            end
            x, τx, gx = xn, τn, gn
            h = xn < absorber_core ? absorber_step_core : xn*(absorber_step_wing - 1)
        end
    end
    [exp(-(depth[k] + mean_depth[k]))*max(0.0, 1 - removed[k]) for k in 1:n]
end

# The logo of Radix: the Chandra observatory points at an accreting black hole and records the Fe K lines of the gas (the spectrum that `spectrum.jl`
# computes with Radix, with noise).
#
#     julia docs/logo/logo.jl
#
# writes `logo.svg` and `logo-dark.svg` to docs/src/assets.
using Printf

const assets = joinpath(@__DIR__, "..", "src", "assets")
const plot_x = (22.0, 222.0)
const baseline = 166.0                     # the y of zero counts
const tallest = 140.0                      # the height of the strongest line (three quarters of the logo or more)

function read_spectrum(path)
    rows = [parse.(Float64, split(line, ',')) for line in readlines(path)[2:end]]
    first.(rows), getindex.(rows, 2), getindex.(rows, 3)
end

xof(e) = plot_x[1] + (e - energy_range[1])/(energy_range[2] - energy_range[1])*(plot_x[2] - plot_x[1])
yof(c) = baseline - c*tallest/peak
f(x) = @sprintf("%.1f", x)

function logo(; ink, spectrum_ink, model, ring, glow, hot, body, panel, panel_line, dark_core, jet, energy, smooth, observed)
    w = (plot_x[2] - plot_x[1])/(length(energy) - 1)
    # the observation: a histogram with a step in each bin
    steps = "M" * f(xof(energy[1]) - w/2) * "," * f(yof(observed[1]))
    for i in eachindex(energy)
        steps *= " H" * f(xof(energy[i]) + w/2)
        i < length(energy) && (steps *= " V" * f(yof(observed[i + 1])))
    end
    ticks = join((@sprintf("<path d=\"M%s,168 L%s,173\" stroke=\"%s\" stroke-opacity=\"0.6\" stroke-width=\"2\" stroke-linecap=\"round\"/>", f(xof(e)), f(xof(e)), ink) for e in (6404.0, 6700.0, 6970.0) if energy_range[1] < e < energy_range[2]))
    axis = @sprintf("<path d=\"M%s,168 L%s,168\" stroke=\"%s\" stroke-opacity=\"0.6\" stroke-width=\"2\" stroke-linecap=\"round\"/>", f(plot_x[1] - 4), f(plot_x[2] + 4), ink)

    # the black hole (upper right): the disk behind it, the shadow, the ring of the light that it bends and the disk in front
    bx, by = 189.0, 46.0
    hole = """
<g transform="translate($bx $by) scale(1.4) translate($(-bx) $(-by)) rotate(-9 $bx $by)">
<ellipse cx="$bx" cy="$by" rx="30" ry="7.5" fill="none" stroke="url(#disk)" stroke-width="6"/>
<path d="M$(bx-16),$(by-1) A16,16 0 0 1 $(bx+16),$(by-1)" fill="none" stroke="$ring" stroke-width="3" stroke-linecap="round" stroke-opacity="0.85"/>
<circle cx="$bx" cy="$by" r="11" fill="$dark_core" stroke="$ring" stroke-width="1.6"/>
<path d="M$(bx-30),$by A30,7.5 0 0 0 $(bx+30),$by" fill="none" stroke="url(#disk)" stroke-width="6"/>
</g>"""

    # Chandra (upper left): the barrel with the mirrors at the front, the science instruments at the back and two solar arrays
    tx, ty, angle = 57.0, 68.0, -9.0
    chandra = """
<g transform="translate($tx $ty) scale(1.4) rotate($angle)">
<path d="M-31,-8 L-17,-8 L-12,-32 L-26,-32 Z" fill="$panel"/>
<path d="M-31,8 L-17,8 L-12,32 L-26,32 Z" fill="$panel"/>
<path d="M-28.5,-8 L-23.5,-32 M-25.5,-8 L-20.5,-32 M-22,-8 L-17,-32 M-18.7,-14.6 L-27,-14.6 M-17.3,-21.2 L-27.7,-21.2 M-15.8,-27.8 L-26.7,-27.8" stroke="$panel_line" stroke-width="0.9" fill="none"/>
<path d="M-28.5,8 L-23.5,32 M-25.5,8 L-20.5,32 M-22,8 L-17,32 M-18.7,14.6 L-27,14.6 M-17.3,21.2 L-27.7,21.2 M-15.8,27.8 L-26.7,27.8" stroke="$panel_line" stroke-width="0.9" fill="none"/>
<rect x="-36" y="-6.5" width="42" height="13" rx="2.5" fill="url(#body)"/>
<rect x="4" y="-8.5" width="32" height="17" rx="3" fill="url(#body)"/>
<rect x="-12" y="-6.5" width="1.6" height="13" fill="$panel_line" fill-opacity="0.5"/>
<ellipse cx="36" cy="0" rx="2.8" ry="8.5" fill="$dark_core" stroke="$hot" stroke-width="1"/>
<rect x="14" y="-8.5" width="1.6" height="17" fill="$panel_line" fill-opacity="0.5"/>
</g>"""
    defs = """<defs>
<linearGradient id="body" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="$(body[1])"/><stop offset="1" stop-color="$(body[2])"/></linearGradient>
<linearGradient id="disk" x1="0" x2="1"><stop offset="0" stop-color="$hot" stop-opacity="0.15"/><stop offset="0.5" stop-color="$glow"/><stop offset="1" stop-color="$hot" stop-opacity="0.15"/></linearGradient>
<linearGradient id="jet" x1="0" y1="1" x2="0" y2="0"><stop offset="0" stop-color="$glow" stop-opacity="0.7"/><stop offset="1" stop-color="$hot" stop-opacity="0"/></linearGradient>
</defs>"""
    """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 240 182"><title>Radix</title>$defs$hole$chandra<path d="$steps" fill="none" stroke="$spectrum_ink" stroke-width="0.9" stroke-linejoin="miter"/>$axis$ticks</svg>
"""
end

energy, smooth, observed = read_spectrum(joinpath(@__DIR__, "fek.csv"))
const energy_range = (energy[1], energy[end])      # eV
const peak = maximum(smooth)                       # the counts of the strongest line of the model
common = (; energy, smooth, observed)
write(joinpath(assets, "logo.svg"), logo(; common..., ink="#24508f", spectrum_ink="#000000", model="#e8833a", ring="#f0a04b", glow="#ffd27a", hot="#e8833a",
    body=("#dfe6f2", "#8e9db8"), panel="#2c4f86", panel_line="#9db8e0", dark_core="#0d1220", jet="#e8833a"))
write(joinpath(assets, "logo-dark.svg"), logo(; common..., ink="#a9c8f5", spectrum_ink="#ffffff", model="#ff9a4d", ring="#ffb35e", glow="#ffe2a8", hot="#ff9a4d",
    body=("#f0f4fb", "#aab8d0"), panel="#4f7bc0", panel_line="#cfe0fb", dark_core="#05070d", jet="#ff9a4d"))

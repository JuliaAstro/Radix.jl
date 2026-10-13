# The energy that Compton scattering exchanges between the radiation and a thermal gas, from first principles and as a surrogate model of it.
#
# XSTAR's `cmpfnc` interpolates linearly in the table `coheat.dat` (101 temperatures × 101 photon energies, a factor 1.23 apart, with the limit 4 Te - E below E = 10⁻⁴). The table
# is not what it should be: compared with the exchange of Klein-Nishina scattering by a relativistic Maxwellian gas (`thermal_exchange`), its entries are good to 0.1% in the typical case but
# wrong by 2% at Te = 10⁻⁴ to 10⁻² and 19% at Te = 0.01 to 1, linear interpolation between them is wrong by up to 76% (the relative error of the value itself is 1 in the 1% of the cells below E = 10⁻⁴ with a
# relativistic electron temperature, where 4 Te - E is no longer the limit), and the cooling of soft spectra by 0.2% (4% at 10⁸ K). `ComptonSurrogate` is a radial basis function surrogate
# (Surrogates.jl) of `thermal_exchange` over the range of the table that matters, which is the default (`compton_surrogate` of the constants); the interpolation of the table is `compton_surrogate = false`.

using Surrogates: RadialBasis, cubicRadial

const compton_Te_range = (1e-9, 10.0)       # electron temperature (electron rest masses) of the surrogate: 6 K to 5×10¹⁰ K
const compton_E_range = (1e-8, 10.0)        # photon energy (electron rest masses): 5 meV to 5 MeV (XSTAR's grid goes from 0.1 eV to 1 MeV)
const compton_surrogate_nodes = (65, 58)    # of the training grid, equally spaced in the logarithms: 6.5 and 6.4 per decade (49 × 44 is 1.7× faster and 8× less accurate; the error falls as N⁻⁴, the time rises as N)
const compton_surrogate_quadrature = 96       # points of each Gauss-Legendre rule of `thermal_exchange` that trains it (64 gives 7×10⁻⁵ of its own error, which the grid of 65 × 58 reaches)
const compton_surrogate_file = joinpath(@__DIR__, "..", "data", "compton_surrogate.dat")

"""
    gauss_legendre(n)

The nodes and weights of the Gauss-Legendre rule of `n` points on [-1, 1].
"""
function gauss_legendre(n::Integer)
    x, w = zeros(n), zeros(n)
    for i in 1:n
        z = cos(π*(i - 1/4)/(n + 1/2))
        dp = 0.0
        for _ in 1:100
            p1, p2 = 1.0, 0.0
            for j in 1:n
                p1, p2 = ((2j - 1)*z*p1 - (j - 1)*p2)/j, p1
            end
            dp = n*(z*p1 - p2)/(z^2 - 1)
            dz = p1/dp
            z -= dz
            abs(dz) < 1e-15 && break
        end
        x[i] = z
        w[i] = 2/((1 - z^2)*dp^2)
    end
    x, w
end

const exchange_rule = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
const exchange_tmax = 30.0                  # kinetic energies of the electrons are integrated up to this many Te (e⁻³⁰)

"""
    thermal_exchange(Te, E; n=64)

The fractional energy that a photon of the energy `E` gains by Compton scattering in a gas of free electrons of the temperature `Te` (both in electron rest masses `mₑc²`), per Thomson cross section:
`⟨ε' - ε⟩/ε` of the scattering rate (Klein-Nishina, so the rate is that of the cross section of the photon in the rest frame of each electron, times `1 - β cos θ`), with the electrons of a Maxwell-Jüttner
distribution. It is XSTAR's `cmpfnc` (the table of `coheat.dat`), `4 Te - E` for `Te, E ≪ 1`. The integral over the Lorentz factor, the angle of the photon to the electron and the angle of the scattering
in the rest frame are done with Gauss-Legendre rules of `n` points (the average over the azimuth is exact: the scattered energy is linear in its cosine); `n = 64` is converged to 10⁻⁴, 96 to 10⁻⁶ of that; the surrogate is trained with 96.
"""
function thermal_exchange(Te::Real, E::Real; n::Integer=64)
    x, w = get!(() -> gauss_legendre(n), exchange_rule, n)
    tmax = exchange_tmax*Te
    num = den = 0.0
    for a in 1:n
        t = tmax*(x[a] + 1)/2                              # the kinetic energy of the electron (mₑc²)
        γ = 1 + t
        β = sqrt(t*(2 + t))/γ
        wγ = γ*sqrt(t*(2 + t))*exp(-t/Te)*tmax/2*w[a]       # the Maxwell-Jüttner weight γ² β e^{-γ/Te} dγ (the constant cancels)
        den += wγ
        for b in 1:n
            μ = x[b]
            flux = 1 - β*μ
            εr = γ*E*flux                                  # the energy of the photon in the rest frame of the electron
            μr = (μ - β)/flux
            for d in 1:n
                cs = x[d]
                r = 1/(1 + εr*(1 - cs))                     # ε'/ε in the rest frame (Compton)
                kn = r^2*(r + 1/r - (1 - cs^2))*π            # Klein-Nishina, in r_e², with the 2π of the azimuth: (r_e²/2) r² (r + 1/r - sin²θ) 2π
                εp = γ*εr*r*(1 + β*μr*cs)                   # back to the frame of the gas, averaged over the azimuth
                s = kn*w[d]*flux*w[b]/2*wγ
                num += s*(εp - E)
            end
        end
    end
    num/den/(8π/3)/E
end

# the surrogate: Es/(Te + E) (bounded, smooth, with no zero where Es changes sign) in the logarithms of Te and E
function train_compton_surrogate(; nodes=compton_surrogate_nodes, n=compton_surrogate_quadrature)
    lt, le = compton_surrogate_axes(nodes)
    values = zeros(nodes...)
    thermal_exchange(1e-3, 1e-3; n)                          # (the rule, before the threads)
    Threads.@threads for k in 1:length(values)
        i, j = Tuple(CartesianIndices(values)[k])
        Te, E = exp(lt[i]), exp(le[j])
        values[i, j] = thermal_exchange(Te, E; n)/(Te + E)
    end
    points = [(a, b) for a in lt, b in le]
    RadialBasis(vec(points), vec(values), (first(lt), first(le)), (last(lt), last(le)); rad=cubicRadial())
end

# the training grid of a surrogate: the nodes of Te and E, equally spaced in the logarithms
compton_surrogate_axes(nodes) = (range(log.(compton_Te_range)...; length=nodes[1]), range(log.(compton_E_range)...; length=nodes[2]))

"""
    save_compton_surrogate(path=compton_surrogate_file; nodes=compton_surrogate_nodes, n=compton_surrogate_quadrature)

Trains the surrogate of `thermal_exchange` and writes what it needs to evaluate it to the text file `path`, as `coheat.dat` is a table: the grid (the number of nodes of Te and E and the two ranges), the
values of `Es/(Te + E)` at the nodes and the coefficients of the radial basis functions (with the other numbers of the model), in 17 digits. `compton_surrogate()` reads it, instead of training the model at each start of Julia.
Run it again if the model changes (`thermal_exchange`, the nodes or the ranges): the tests compare the file with a new training.
"""
function save_compton_surrogate(path::AbstractString=compton_surrogate_file; nodes=compton_surrogate_nodes, n=compton_surrogate_quadrature)
    s = train_compton_surrogate(; nodes, n)
    open(path, "w") do io
        println(io, "# Radix: the radial basis function surrogate (Surrogates.jl, cubic) of the energy exchanged in Compton scattering (`thermal_exchange`),")
        println(io, "# as Es/(Te + E) on the grid of the logarithms of Te and E, with its coefficients; made by `Radix.save_compton_surrogate`")
        println(io, nodes[1], " ", nodes[2], " ", compton_Te_range[1], " ", compton_Te_range[2], " ", compton_E_range[1], " ", compton_E_range[2])
        println(io, s.dim_poly, " ", s.scale_factor, " ", length(s.coeff))
        foreach(v -> println(io, repr(v)), s.y)
        foreach(v -> println(io, repr(v)), vec(s.coeff))
    end
    path
end

"""
    load_compton_surrogate(path=compton_surrogate_file)

The surrogate of `thermal_exchange` that `save_compton_surrogate` wrote to `path`, without training it. (If this version of Surrogates.jl does not take the coefficients, the system for them is solved again from the values, which takes 0.2 s.)
"""
function load_compton_surrogate(path::AbstractString=compton_surrogate_file)
    lines = [l for l in eachline(path) if !startswith(l, "#") && !isempty(strip(l))]
    h = split(lines[1])
    nodes = (parse(Int, h[1]), parse(Int, h[2]))
    ranges = parse.(Float64, h[3:6])
    (ranges[1], ranges[2]) == compton_Te_range && (ranges[3], ranges[4]) == compton_E_range ||
        error("the surrogate in $path is for other ranges of Te and E than those of Radix")
    p = split(lines[2])
    dim_poly, scale_factor, ncoeff = parse(Int, p[1]), parse(Float64, p[2]), parse(Int, p[3])
    nvalues = prod(nodes)
    values = parse.(Float64, lines[3:2 + nvalues])
    coeff = parse.(Float64, lines[3 + nvalues:2 + nvalues + ncoeff])
    lt, le = compton_surrogate_axes(nodes)
    points = vec([(a, b) for a in lt, b in le])
    lb, ub = (first(lt), first(le)), (last(lt), last(le))
    rad = cubicRadial()
    try
        RadialBasis(rad.phi, dim_poly, points, values, lb, ub, transpose(coeff), scale_factor, false, 0.0)
    catch
        RadialBasis(points, values, lb, ub; rad=cubicRadial())
    end
end

const compton_surrogate_cache = Ref{Any}(nothing)
const compton_surrogate_lock = ReentrantLock()

"""
    compton_surrogate()

The radial basis function surrogate of `thermal_exchange` of the Compton scattering, read at the first call from `compton_surrogate_file` (`load_compton_surrogate`; 3770 points: a few milliseconds). If the file is not there it is trained (`train_compton_surrogate`: 1-10 s)
and a warning says so.
"""
function compton_surrogate()
    s = compton_surrogate_cache[]
    s === nothing || return s
    lock(compton_surrogate_lock) do
        if compton_surrogate_cache[] === nothing
            compton_surrogate_cache[] = if isfile(compton_surrogate_file)
                load_compton_surrogate()
            else
                @warn "No $compton_surrogate_file: training the surrogate of the Compton scattering"
                train_compton_surrogate()
            end
        end
        compton_surrogate_cache[]
    end
end

"""
    surrogate_exchange(E, Te)

`thermal_exchange` from the surrogate, or `nothing` outside of its range of `Te` and `E` (above 10 `mₑc²` of either; the table is used there). Below the range the limit `4 Te - E` is added to the value at its edge.
"""
function surrogate_exchange(E::Real, Te::Real)
    (E > compton_E_range[2] || Te > compton_Te_range[2]) && return nothing
    Tc, Ec = max(Te, compton_Te_range[1]), max(E, compton_E_range[1])
    compton_surrogate()((log(Tc), log(Ec)))*(Tc + Ec) + 4*(Te - Tc) - (E - Ec)
end

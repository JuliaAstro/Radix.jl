# The radiation field on XSTAR's energy grid, and the continuum opacity and
# emissivity arrays that the photoionization integrals add to.

"""
    Radiation(E, F)

Photon energy grid `E` (eV, geometrically spaced; XSTAR's `epi`) and the incident
spectrum `F` on it (XSTAR's `bremsa`). The photoionization rates integrate cross
sections over it. `Fint` (XSTAR's `bremsint`) is the flux integrated from each bin to
the top of the grid, in erg.
"""
struct Radiation{R<:AbstractFloat}
    E::Vector{R}
    F::Vector{R}
    Fint::Vector{R}
    function Radiation(E::AbstractVector{R}, F::AbstractVector{R}) where {R<:AbstractFloat}
        length(E) == length(F) || throw(DimensionMismatch("E and F must have the same length"))
        n = length(E)
        Fint = zeros(R, n)
        ergsev = constants().ergsev_bremsint
        for k in n - 1:-1:1
            Fint[k] = Fint[k + 1] + (F[k] + F[k + 1])*(E[k + 1] - E[k])/2*ergsev
        end
        new{R}(collect(E), collect(F), Fint)
    end
end

"""
    Opacity([R=Float64,] n)

Continuum opacity arrays on a grid of `n` energies that the photoionization rates
add to when given as `opacity=`: `total` (XSTAR `opakc`, cm⁻¹), `continuum`
(`opakcont`, lines excluded) and `emissivity` (`rccemis`, 2 × n, recombination
continuum emissivities in and out). `R` is the element type: use the type of the temperature
when differentiating the rates with ForwardDiff.
"""
struct Opacity{R<:Real}
    total::Vector{R}
    continuum::Vector{R}
    emissivity::Matrix{R}
end
Opacity(R::Type{<:Real}, n::Integer) = Opacity(zeros(R, n), zeros(R, n), zeros(R, 2, n))
Opacity(n::Integer) = Opacity(Float64, n)

# grid constants of nbinc/huntf
const grid_guard_fraction = 50       # the top 1/50 of the grid (at least 2 bins) is not used
const grid_min_bins = 2
const grid_tiny = 1e-34

"""
    nbin(rad, e)

Index of the grid energy nearest to `e` (XSTAR's `nbinc`): the grid is assumed
equally spaced in log E.
"""
function nbin(rad::Radiation, e)
    n = length(rad.E) - max(grid_min_bins, length(rad.E) ÷ grid_guard_fraction)
    xx = rad.E
    xtmp = max(e, xx[2])
    (e < grid_tiny || xx[1] <= grid_tiny || xx[n] <= grid_tiny) && return 1
    jlo = trunc(Int, (n - 1)*log(xtmp/xx[1])/log(xx[n]/xx[1])) + 1
    if jlo < n
        abs(log(e/(grid_tiny + xx[jlo + 1]))) < abs(log(e/(grid_tiny + xx[jlo]))) && (jlo += 1)
    end
    clamp(jlo, 1, n)
end

"""
    xstar_energy_grid(ncn2=9999)

XSTAR's photon energy grid (its `ener`), in eV: `ncn2 - max(2, ncn2/50)` points spaced geometrically from 0.1 to
4×10⁵ eV, then `max(2, ncn2/50)` more from there to 10⁶ eV.
"""
function xstar_energy_grid(ncn2=9999)
    numcon2 = max(grid_min_bins, ncn2 ÷ grid_guard_fraction)
    numcon3 = ncn2 - numcon2
    E = zeros(ncn2)
    E[1] = 0.1
    dele = (4e5/0.1)^(1/(numcon3 - 1))
    for l in 2:numcon3
        E[l] = E[l - 1]*dele
    end
    dele = (1e6/4e5)^(1/(numcon2 - 1))
    for l in numcon3 + 1:ncn2
        E[l] = E[l - 1]*dele
    end
    E
end

"""
    NO_RADIATION

A radiation field without photons on XSTAR's energy grid: the default of the photoionization rates, for which the
photoionization rates are 0 and the recombination rates are those of the gas alone.
"""
const NO_RADIATION = let E = xstar_energy_grid()
    Radiation(E, zeros(length(E)))
end

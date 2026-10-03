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
        for k in n - 1:-1:1
            Fint[k] = Fint[k + 1] + (F[k] + F[k + 1])*(E[k + 1] - E[k])/2*ergsev_bremsint
        end
        new{R}(collect(E), collect(F), Fint)
    end
end

"""
    Opacity(n)

Continuum opacity arrays on a grid of `n` energies that the photoionization rates
add to when given as `opacity=`: `total` (XSTAR `opakc`, cm⁻¹), `continuum`
(`opakcont`, lines excluded) and `emissivity` (`rccemis`, 2 × n, recombination
continuum emissivities in and out).
"""
struct Opacity{R<:AbstractFloat}
    total::Vector{R}
    continuum::Vector{R}
    emissivity::Matrix{R}
end
Opacity(n::Integer) = Opacity(zeros(n), zeros(n), zeros(2, n))

const ergsev_bremsint = 1.602197e-12     # erg per eV in bremsmap's integral

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

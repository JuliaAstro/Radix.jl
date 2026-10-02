# Local gas state of a cell, in XSTAR's units
#   T:    temperature in units of 10⁴ K
#   nₕ:   neutral hydrogen density (cm⁻³), the XSTAR `xh0`
#   nₑ:   electron density (cm⁻³), the XSTAR `xnx`
#   ntot: total hydrogen density (cm⁻³), the XSTAR `xpx`
#
# Static atomic data (level energies and weights, element masses) live in tables
# and the radiation field is a separate object, not part of Cell.

struct Cell{R}
    T::R
    nₕ::R
    nₑ::R
    ntot::R
end

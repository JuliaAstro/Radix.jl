# Physical conditions of a cell, in XSTAR's units
#   T:  temperature in units of 10⁴ K
#   nₕ: neutral hydrogen density (cm⁻³), the XSTAR `xh0`
#   nₑ: electron density (cm⁻³), the XSTAR `xnx`

struct Cell{R}
    T::R
    nₕ::R
    nₑ::R
end

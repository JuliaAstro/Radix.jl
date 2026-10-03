# Physical constants and XSTAR conventions used by more than one file. Numbers that
# belong to a single fit or rate are named in the file that uses them.

const mc² = 5.11e5                    # electron rest energy (eV)
const ergsev = 1.602176634e-12        # erg per eV (xstarlib constants module)
const Ry_eV = 13.605692               # Rydberg energy (eV), as used in ucalc
const Ry_eV_single = Float64(Float32(Ry_eV))   # 13.605692 as the single-precision literal of ucalc
const hc_eVÅ = 12398.4016             # h c (eV Å), as used in ucalc

# XSTAR measures temperature in units of 10⁴ K
const T_unit = 1e4                    # K per unit of temperature
const kT_eV = 0.861707                # k_B × 10⁴ K (eV)
const T_floor_coeff = 2.8777e6        # K Å: the fit temperature of the collision
                                      # rates is at least T_floor_coeff/λ (ΔE/kT ≤ 50)

# Maxwellian collision rates: C = collision_rate_coeff Υ / (√T g) in cm³ s⁻¹
const collision_rate_coeff = 8.626e-8

const cx_unit = 1e-9                  # charge-exchange rates are fitted in 10⁻⁹ cm³ s⁻¹
const tiny = 1e-24                    # guard against division by zero and zero wavelengths

# exp() with the argument clamped, as in XSTAR's expo()
const expo_limit = 60.0
expo(x) = exp(clamp(x, -expo_limit, expo_limit))

const kB_cgs = 1.380649e-16           # Boltzmann constant (erg K⁻¹)
const fourpi_xstar = 12.56            # XSTAR's value of 4π in the spectrum integrals
const Mb = 1e-18                      # cm² per Mb (cross sections are tabulated in Mb)

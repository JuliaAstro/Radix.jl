# The physical constants the rates use, from PhysicalConstants.jl (CODATA 2022 by default), and
# the rounded values of XSTAR's `ucalc` that the reference tests need.

using PhysicalConstants
using Unitful
import PhysicalConstants.CODATA2022

"""
    Constants

The physical constants and the derived combinations the rates use, in the units of the rates (cgs, eV,
Å, K). `Constants()` (or `Constants(PhysicalConstants.CODATA2018)`, ...) evaluates them from a CODATA
set of PhysicalConstants.jl, `Constants(K; kT_eV=...)` copies `K` with some fields changed, `ucalc_constants()` has the rounded values written in XSTAR's source.

`ucalc` rounds the same constant differently in different routines (its k is 8.617e-5, 0.861707e-4
and 1.38066e-16 in places), so some physical constants have more than one field, one per use; they
are equal for CODATA. Fields:

- `ergsev` erg per eV, with `ergsev_bremsint` (the spectrum integral) and `ergsev_decay` (type 71) variants.
- `kB_cgs` (erg K⁻¹), `kB_szirc`, `kB_eV` (eV K⁻¹) and `kT_eV` (k × 10⁴ K, in eV): Boltzmann's constant.
- `Ry_eV` (hcR∞), `Ry_eV_single` (`ucalc`'s single-precision literal), `Ry_erg` and `Ry_erg_szirc`.
- `hc_eVÅ` (eV Å) with `hc_eVÅ_single`, `hc_over_k` (Å K) and `T_floor_coeff` (K Å; the collision rates
  use a temperature of at least `T_floor_coeff/λ`, i.e. ΔE/50k).
- `Ry_K`, `Ry_K_erc`, `Ry_K_cb` (the Rydberg over k, in K), `eV_K` (eV/k in K), `inv_kB` (K erg⁻¹).
- `invcm_per_eV` (cm⁻¹ of 1 eV), `Rinf_invcm` (cm⁻¹).
- `fourpi`, `eightpi` and `pi_pexs`: 4π and 8π (the normalization of the spectrum in the photoionization
  integrals) and π (the resonance profiles of type 85). These are mathematical constants, not CODATA,
  but `ucalc` rounds them (12.56, 25.3 and 3.14159).
- Coefficients that are combinations of constants in XSTAR's source: `bb_coeff` (1/(h³c²), eV⁻³ s⁻¹ cm⁻²),
  `fo_saha` (`saha_coeff` × 8π with T in 10⁴ K), `saha_coeff_cb` (the same factor in the ionization
  of type 57), `A_to_f` (1/(8π² r_e c) in Å² s⁻¹), `milne_coeff` (4π/((2π m_e)^{3/2} c²) × Mb, cgs),
  `gordon_A` (the hydrogenic dipole rate, (2π/3) α³ c R∞ s⁻¹ with the reduced mass of hydrogen), `ps_alfa` (m_e/2k), `ps_pd` (the
  Debye length coefficient (k/(4π α ħ c))^{1/2}), `Ry_K_sz` (Ry/k in the Simpson-Zhang excitation)
  and `sz_rate_coeff` (the Maxwellian rate coefficient with T in K).
- `proton_mass` in electron masses, `light_speed` (cm s⁻¹), `line_xsec` (π e²/(m_e c), cm² Hz),
  `thermal_speed` (cm s⁻¹ of an atom of 1 amu at 10⁴ K),
  `collision_rate_coeff` (Maxwellian Υ rate coefficient for T in 10⁴ K, cm³ s⁻¹) and `saha_coeff` (cm³ K^{3/2}).
"""
struct Constants
    ergsev::Float64
    ergsev_bremsint::Float64
    ergsev_decay::Float64
    kB_cgs::Float64
    kB_szirc::Float64
    kB_eV::Float64
    kT_eV::Float64
    Ry_eV::Float64
    Ry_eV_single::Float64
    Ry_erg::Float64
    Ry_erg_szirc::Float64
    hc_eVÅ::Float64
    hc_eVÅ_single::Float64
    hc_over_k::Float64
    T_floor_coeff::Float64
    Ry_K::Float64
    Ry_K_erc::Float64
    Ry_K_cb::Float64
    eV_K::Float64
    inv_kB::Float64
    invcm_per_eV::Float64
    Rinf_invcm::Float64
    fourpi::Float64
    eightpi::Float64
    pi_pexs::Float64
    bb_coeff::Float64
    fo_saha::Float64
    saha_coeff_cb::Float64
    A_to_f::Float64
    milne_coeff::Float64
    gordon_A::Float64
    ps_alfa::Float64
    ps_pd::Float64
    Ry_K_sz::Float64
    sz_rate_coeff::Float64
    proton_mass::Float64
    light_speed::Float64
    line_xsec::Float64
    thermal_speed::Float64
    collision_rate_coeff::Float64
    saha_coeff::Float64
end

const T_floor_dE_over_kT = 50.0       # the fit temperature of the collision rates is at least ΔE/(50 k)
const T_unit = 1e4                    # XSTAR measures temperature in units of 10⁴ K (K)

# a Constants from a NamedTuple of all the fields, or from another Constants with some changed
Constants(values::NamedTuple) = Constants((values[name] for name in fieldnames(Constants))...)
Constants(K::Constants; changes...) =
    Constants(merge(NamedTuple{fieldnames(Constants)}(getfield.(Ref(K), fieldnames(Constants))), values(changes)))

function Constants(codata::Module=CODATA2022)
    k, h, c = codata.BoltzmannConstant, codata.PlanckConstant, codata.SpeedOfLightInVacuum
    mₑ, mₚ, mᵤ = codata.ElectronMass, codata.ProtonMass, codata.AtomicMassConstant
    R∞, a₀, α, ħ = codata.RydbergConstant, codata.BohrRadius, codata.FineStructureConstant,
        codata.ReducedPlanckConstant
    value(unit, x) = Float64(ustrip(unit, x))
    eV = 1u"eV"
    ergsev = value(u"erg", eV)
    kB = value(u"erg/K", k)
    Ry = value(u"eV", h*c*R∞)
    Ry_erg = value(u"erg", h*c*R∞)
    hc = value(u"eV*angstrom", h*c)
    hc_k = value(u"angstrom*K", h*c/k)
    rₑ = α^2*a₀                                       # classical electron radius
    saha = value(u"cm^3*K^(3/2)", (h^2/(2π*mₑ*k))^1.5/2)
    collision = value(u"cm^3/s*K^(1/2)", sqrt(2π)*ħ^2/(mₑ^1.5*k^0.5))
    Constants((;
        ergsev, ergsev_bremsint=ergsev, ergsev_decay=ergsev,
        kB_cgs=kB, kB_szirc=kB, kB_eV=value(u"eV/K", k), kT_eV=value(u"eV/K", k)*T_unit,
        Ry_eV=Ry, Ry_eV_single=Ry, Ry_erg, Ry_erg_szirc=Ry_erg,
        hc_eVÅ=hc, hc_eVÅ_single=hc, hc_over_k=hc_k, T_floor_coeff=hc_k/T_floor_dE_over_kT,
        Ry_K=value(u"K", h*c*R∞/k), Ry_K_erc=value(u"K", h*c*R∞/k), Ry_K_cb=value(u"K", h*c*R∞/k),
        eV_K=value(u"K", eV/k), inv_kB=value(u"K/erg", 1/k),
        invcm_per_eV=value(u"cm^-1", eV/(h*c)), Rinf_invcm=value(u"cm^-1", R∞),
        fourpi=4π, eightpi=8π, pi_pexs=π,
        bb_coeff=value(u"eV^-3*s^-1*cm^-2", 1/(h^3*c^2)),
        fo_saha=saha*T_unit^-1.5*8π, saha_coeff_cb=saha,
        A_to_f=value(u"s/angstrom^2", 1/(8π^2*rₑ*c)),
        milne_coeff=value(u"g^(-3/2)*cm^-2*s^2", 4π/((2π*mₑ)^1.5*c^2))*Mb,
        gordon_A=value(u"Hz", 2π/3*α^3*c*R∞/(1 + mₑ/mₚ)), ps_alfa=value(u"K*s^2/cm^2", mₑ/(2k)),
        ps_pd=value(u"cm^(-1/2)*K^(-1/2)", sqrt(k/(4π*α*ħ*c))),
        Ry_K_sz=value(u"K", h*c*R∞/k), sz_rate_coeff=collision,
        proton_mass=value(Unitful.NoUnits, mₚ/mₑ), light_speed=value(u"cm/s", c),
        line_xsec=value(u"cm^2/s", π*rₑ*c),
        thermal_speed=value(u"cm/s", sqrt(2*k*T_unit*u"K"/mᵤ)),
        collision_rate_coeff=collision/sqrt(T_unit), saha_coeff=saha))
end

"""
    ucalc_constants()

The values XSTAR's `ucalc` uses, rounded and not always consistent with each other. Reproduces its
numbers to about 1e-6; see `with_constants`.
"""
ucalc_constants() = Constants((;
    ergsev=1.602176634e-12, ergsev_bremsint=1.602197e-12, ergsev_decay=Float64(1.602197f-12),
    kB_cgs=1.380649e-16, kB_szirc=Float64(1.38066f-16), kB_eV=8.617e-5, kT_eV=0.861707,
    Ry_eV=13.605692, Ry_eV_single=Float64(13.605692f0), Ry_erg=2.17896e-11,
    Ry_erg_szirc=Float64(2.179874f-11),
    hc_eVÅ=12398.4016, hc_eVÅ_single=Float64(12398.4016f0), hc_over_k=Float64(1.43817f8),
    T_floor_coeff=2.8777e6,
    Ry_K=157888.0, Ry_K_erc=157803.0, Ry_K_cb=Float64(13.605692f0*1.6021f-19/1.3805f-23),
    eV_K=Float64(1.16058f4), inv_kB=7.2438e15, invcm_per_eV=8065.48, Rinf_invcm=109737.0,
    fourpi=12.56, eightpi=25.3, pi_pexs=3.14159,
    bb_coeff=1.571e22, fo_saha=5.216e-21, saha_coeff_cb=Float64(2.0779f-16), A_to_f=1e-16/0.667274,
    milne_coeff=Float64(0.79788f0*40.4153f0), gordon_A=2.6761e9, ps_alfa=3.297e-12, ps_pd=6.90,
    Ry_K_sz=1.578203e5, sz_rate_coeff=8.63e-6,
    proton_mass=1800.0, light_speed=3e10, line_xsec=0.02655, thermal_speed=1.29e6,
    collision_rate_coeff=8.626e-8, saha_coeff=2.07e-16))

const CURRENT_CONSTANTS = Ref(Constants())

"""
    constants()

The set of constants the rates are using (CODATA 2022 unless changed).
"""
constants() = CURRENT_CONSTANTS[]

"""
    set_constants!(K::Constants)

Make `K` the constants of all the rates, e.g. `set_constants!(Constants(PhysicalConstants.CODATA2018))`.
Returns `K`. This is global state: use `with_constants` for a temporary change.
"""
set_constants!(K::Constants) = (CURRENT_CONSTANTS[] = K)

"""
    with_constants(f, K::Constants)

Call `f()` with `K` as the constants of the rates, restoring the previous set afterwards.
"""
function with_constants(f, K::Constants)
    old = CURRENT_CONSTANTS[]
    CURRENT_CONSTANTS[] = K
    try
        f()
    finally
        CURRENT_CONSTANTS[] = old
    end
end

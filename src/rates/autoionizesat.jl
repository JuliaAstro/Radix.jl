# Autoionization rates for satellite levels: r1 = Aa(k,i) (s−1); r2 = E(k) (eV
# above ionization limit); r3= (2J + 1); i1= (2S + 1); i2= L; i3= k (level);
# i4= i (continuum level); i5= Z; i6= ionN ; s1= level configuration

# XSTAR data type: 72

const AutoionizeSatDesc = "Autoinization rates (in s^-1) for satellite lvls"

struct AutoionizeSat{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    spin_mult::I        # 2S+1
    L::I                # orbital angular momentum
    level::I            # level index
    parent::Parent{I}           # parent ion (0: not stored) and continuum level
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    A_auto::R           # autoionization rate (s⁻¹)
    E::R                # energy above ionization limit (eV)
    g::R                # statistical weight 2J+1
end

function AutoionizeSat(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AutoionizeSat(Int8(rate), label, ivec[1:3]..., Parent(zero(eltype(ivec)), ivec[4]),
        ivec[5:6]..., rvec...)
end

const sat_capture_coeff = 3.3e-11      # cm³ s⁻¹ at kT = 13.6 eV
const sat_capture_scale = 1e13         # the stored rate is divided by this

# calt72: dielectronic-capture style rate (s⁻¹) from the stored autoionization rate A (s⁻¹), the energy E above the limit
# (eV) and the weight g (1 if the record has none), at the temperature T (10⁴ K)
function satellite_capture(A, E, g, T)
    K = constants()
    kT = K.kT_eV*T
    sat_capture_coeff*(K.Ry_eV_coarse/kT)^1.5*exp(-E/kT)*(A/sat_capture_scale)*g
end

"""
    rate(coef::AutoionizeSat, cell; levels, nlev, index=false)

Autoionization of a satellite level (XSTAR ucalc type 72) from the rate `3.3e-11 (13.6 eV/kT)^{3/2} e^{-E/kT} (A/10¹³) g`
(`calc_72`). As `ucalc` has it, `irate` is that rate times the electron density and `frate` is `irate nₑ` times the
Saha factor of the ground to the continuum level weights and `e^{E/T}` with the energy in eV divided by the
temperature in K (which is ≈ 1; the units look wrong in the source and are reproduced). `init` is the level and
`final` the continuum level of the record. `levels` needs the ground (1) and the continuum (`nlev`) levels.
"""
function rate(coef::AutoionizeSat, cell::Cell; levels, nlev, index=false, verbose=false)
    K = constants()
    init, final = Int(coef.level), Int(coef.parent.level)
    index && return (; init, final, frate=0., irate=0.)
    ground = get(levels, (coef.ion, 1), nothing)
    cont = get(levels, (coef.ion, nlev), nothing)
    (ground === nothing || cont === nothing) && return (; init, final, frate=0., irate=0.)
    capture = satellite_capture(Float64(coef.A_auto), Float64(coef.E), Float64(coef.g), cell.T)
    irate = capture*cell.nₑ
    rinf = K.saha_ci*Float64(ground.g)/Float64(cont.g)/cell.T/sqrt(cell.T)
    frate = irate*rinf*cell.nₑ*expo(Float64(coef.E)/(cell.T*T_unit))
    (; init, final, frate, irate)
end

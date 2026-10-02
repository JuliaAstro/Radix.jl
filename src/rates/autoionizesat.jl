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
    continuum_level::I  # continuum level
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    A_auto::R           # autoionization rate (s⁻¹)
    E::R                # energy above ionization limit (eV)
    g::R                # statistical weight 2J+1
end

function AutoionizeSat(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AutoionizeSat(Int8(rate), label, ivec..., rvec...)
end

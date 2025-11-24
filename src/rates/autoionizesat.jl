# Autoionization rates for satellite levels: r1 = Aa(k,i) (s−1); r2 = E(k) (eV
# above ionization limit); r3= (2J + 1); i1= (2S + 1); i2= L; i3= k (level);
# i4= i (continuum level); i5= Z; i6= ionN ; s1= level configuration

# XSTAR data type: 72

const AutoionizeSatDesc = "Autoinization rates (in s^-1) for satellite lvls"

struct AutoionizeSat{I,R} <: AbstractRate
    rate::Int8
    label::String
    S2::I
    L::I
    k::I
    i::I
    Z::I
    N::I
    A::R
    E::R
    J2::R
end

function AutoionizeSat(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AutoionizeSat(Int8(rate), label, ivec..., rvec...)
end

function rate(coef::AutoionizeSat, cell::Cell)
    
end

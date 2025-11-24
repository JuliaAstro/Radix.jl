# Autoionization rates for Fe XXIV satellites [8]: r1 = Aa(k,i) (s−1);
# r2 = E(k) (eV above ionization limit); i1= ionN , i2= kN ; i3= ionN−1;
# i4= iN−1; i5= ionN

# XSTAR data type: 75

const AutoionizeFe25SatDesc = "autoionization data for Fe XXiV satellites"

struct AutoionizeFe25Sat{I,R} <: AbstractRate
    rate::Int8
    label::String
    N::I
    k::I
    Nm::I
    i::I
    N_::I
    A::R
    E::R
end

function AutoionizeFe25Sat(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AutoionizeFe25Sat(Int8(rate), label, ivec..., rvec...)    
end

function rate(coef::AutoionizeFe25Sat, cell::Cell)
    
end

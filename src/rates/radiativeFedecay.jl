# Decay rates for Fe UTA [24]: r1 = λ(Å); r2 = E(k) (eV); r3= g f (i,k);
# r4 = Ar(k,i) (s−1); r5 = Aa(k,i) (s−1); i1= i (lower level);
# i2= k (upper level); i4= ionN

# XSTAR data type: 82

const RadiativeFeDecayDesc = "Fe UTA rad rates"

struct RadiativeFeDecay{I,R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    N::I
    λ::R
    E::R
    g::R
    Ar::R
    Aa::R
end

function RadiativeFeDecay(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    RadiativeFeDecay(Int8(rate), label, ivec..., rvec...)
end

function rate(coef::RadiativeFeDecay)
    
end

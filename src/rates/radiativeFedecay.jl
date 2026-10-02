# Decay rates for Fe UTA [24]: r1 = λ(Å); r2 = E(k) (eV); r3= g f (i,k);
# r4 = Ar(k,i) (s−1); r5 = Aa(k,i) (s−1); i1= i (lower level);
# i2= k (upper level); i4= ionN

# XSTAR data type: 82

const RadiativeFeDecayDesc = "Fe UTA rad rates"

struct RadiativeFeDecay{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    ion::I     # ion index (XSTAR ionN)
    λ::R       # wavelength (Å)
    E::R       # eV
    gf::R      # weighted oscillator strength
    A_rad::R   # s⁻¹
    A_auto::R  # s⁻¹
end

function RadiativeFeDecay(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    RadiativeFeDecay(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3], rvec...)
end

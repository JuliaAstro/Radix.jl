# Auger and radiative widths of kN -th K-vacancy level: r1 = E(kN ) (eV,
# relative to E(∞)); r2 = Aa(kN ) (s−1); r3= Aa(kN ,iN−1) (s−1);
# r4 = Ar(kN ) (s−1); i1= iN−1; i2= kN ; i3= Z; i4= ionN−1; i5= ionN

# XSTAR data type: 86

const IronKAugerDesc = "Iron K Auger data from Patrick"

struct IronKAuger{I,R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    Nm::I
    N::I
    E::R
    A::Tuple
end

function IronKAuger(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    IronKAuger(Int8(rate), label, ivec..., rvec[1], Tuple(rvec[2:end]))
end

function rate(coef::IronKAuger)

end

# Auger and radiative widths of kN -th K-vacancy level: r1 = E(kN ) (eV,
# relative to E(∞)); r2 = Aa(kN ) (s−1); r3= Aa(kN ,iN−1) (s−1);
# r4 = Ar(kN ) (s−1); i1= iN−1; i2= kN ; i3= Z; i4= ionN−1; i5= ionN

# XSTAR data type: 86

const IronKAugerDesc = "Iron K Auger data from Patrick"

struct IronKAuger{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    parent::Parent{I}           # parent ion and level
    level::I             # level index
    Z::I                 # atomic number
    ion::I               # ion index (XSTAR ionN)
    E::R                 # energy relative to E∞ (eV)
    A_widths::Vector{R}  # [A_auto(k), A_auto(k,parent), A_rad(k)] (s⁻¹)
end

function IronKAuger(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    IronKAuger(Int8(rate), label, Parent(ivec[4], ivec[1]), ivec[2], ivec[3], ivec[5], rvec[1], Vector(rvec[2:end]))
end

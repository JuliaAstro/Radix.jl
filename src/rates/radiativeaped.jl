# APED line (k−i) radiation rates [20]: r1 = λ(Å); r2 = 0.0; r3= A(k,i) (s−1);
# i1= i (lower level); i2= k (upper level); i3= Z; i4= ionN

# XSTAR data type: 91

const RadiativeAPEDDesc = "aped line wavelengths same as 50"

struct RadiativeAPED{I, R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    N::I
    λ::R
    A::R
end

function RadiativeAPED(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    RadiativeAPED(Int8(rate), label, ivec..., rvec[1], rvec[3])
end

function rate(coef::RadiativeAPED)
    
end

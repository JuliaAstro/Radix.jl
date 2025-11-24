# Line (k−i) radiation rates of N-electron ion: r1 = λ(Å); r2 = g f (i,k);
# r3= A(k,i) (s−1); i1= i (lower level); i2= k (upper level); i3= Z; i4= ionN

# XSPEC data type: 50

const AtomicLine2Desc = "op line rad. rates"

struct AtomicLine2{I, R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    N::I
    λ::R
    gf::R
    A::R
end

function AtomicLine2(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AtomicLine2(Int8(rate), label::String, ivec..., rvec...)
end

function rate(coef::AtomicLine2, cell::Cell)

end

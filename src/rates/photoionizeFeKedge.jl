# Photoionization cross sections for Fe ions obtained by summation of
# resonances near the K edge [44]: r1 = Zeff; r2 = Eth (Ryd); r3= f; r4 = γ;
# r5 = scaling factor; i1= n; i2= L; i3= 2J; i4= Z; i5= kN−1; i6= ionN−1;
# i7= iN ; i8= ionN

# XSTAR data type: 85

const PhotoionizeFeKedgeDesc = "Iron K Pi xsections, spectator Auger summed"

struct PhotoionizeFeKedge{I,R} <: AbstractRate
    rate::Int8
    label::String
    n::I
    L::I
    J2::I
    Z::I
    km::I
    Nm::I
    i::I
    N::I
    Zeff::R
    E::R
    f::R
    γ::R
    c::R
end

function PhotoionizeFeKedge(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    PhotoionizeFeKedge(Int8(rate), label, ivec..., rvec...)
end

function rate(coef::PhotoionizeFeKedge)
    
end

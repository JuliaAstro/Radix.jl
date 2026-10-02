# Photoionization cross sections for Fe ions obtained by summation of
# resonances near the K edge [44]: r1 = Zeff; r2 = Eth (Ryd); r3= f; r4 = γ;
# r5 = scaling factor; i1= n; i2= L; i3= 2J; i4= Z; i5= kN−1; i6= ionN−1;
# i7= iN ; i8= ionN

# XSTAR data type: 85

const PhotoionizeFeKedgeDesc = "Iron K Pi xsections, spectator Auger summed"

struct PhotoionizeFeKedge{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I             # principal quantum number
    L::I             # orbital angular momentum
    twoJ::I          # 2J
    Z::I             # atomic number
    parent::Parent{I}           # parent ion and level
    level::I         # level index
    ion::I           # ion index (XSTAR ionN)
    Zeff::R
    E_th::R          # Ry
    f::R
    γ::R
    scale::R         # scaling factor
end

function PhotoionizeFeKedge(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    PhotoionizeFeKedge(Int8(rate), label, ivec[1:4]..., Parent(ivec[6], ivec[5]), ivec[7:8]..., rvec...)
end

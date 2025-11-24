# Photoionization cross section damped excess of iN -th level of the N-electron
# ion leaving the (N−1)-electron ion in superlevel_[K] kN−1: r1−rjmax = 
# (E( j),σ(E( j)), j= 1, jmax) (Energy in Ryd relative to E(∞), cross section
# in Mb); i1= n; i2= L; i3= 2J; i4= Z; i5= kN−1; i6= iN ; i7= ionN.

# XSTAR data type: 88

const PhotoionizeDampDesc = "Iron inner shell resonance excitation (Patrick)"

struct PhotoionizeDamp{I} <: AbstractRate
    rate::Int8
    label::String
    n::I
    L::I
    J2::I
    Z::I
    kNm::I
    k::I
    i::I
    N::I
    E::Tuple
    σ::Tuple
end

function PhotoionizeDamp(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    iv = length(ivec) < 8 ? (ivec[1:5]..., Int32(0), ivec[6:7]...) : Tuple(ivec)
    PhotoionizeDamp(Int8(rate), label, iv..., Tuple(rvec[1:2:end-1]),
        Tuple(rvec[2:2:end]))
end

function rate(coef::PhotoionizeDamp)

end

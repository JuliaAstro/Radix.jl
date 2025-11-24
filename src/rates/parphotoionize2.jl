# TOPbase partial photoionization cross section (resonance averaged) of iN-th
# level of the N-electron ion leaving the (N−1)-electron ion in the kN−1-th
# level: r1−rj2*max = (E( j),σ(E( j)), j= 1, jmax) (Energy in Ryd relative to
# E(∞), cross section in Mb); i1= n; i2= L; i3= 2J; i4= Z; i5= kN−1;
# i6= ionN−1; i7= iN ; i8= ionN

# XSTAR data type: 53

const ParPhotoIonize2Desc = "op pi xsections"

struct ParPhotoIonize2{I} <: AbstractRate
    rate::Int8
    label::String
    n::I
    L::I
    J2::I
    Z::I
    kNm::I
    ionNm::I
    iN::I
    N::I
    E::Tuple
    σ::Tuple
end

function ParPhotoIonize2(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ParPhotoIonize2(Int8(rate), label, ivec..., Tuple(rvec[1:2:end-1]),
        Tuple(rvec[2:2:end]))
end

function rate(coef::ParPhotoIonize2)

end

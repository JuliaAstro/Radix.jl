# Electron-impact effective collision strengths for the k−i transition of
# N-electron ion: r1−rjmax = (log Te( j), j= 1, jmax) (K); rj(max+1)−rj(2*max)=
# (Υ(Te( j)), j= 1, jmax) (effective collision strength); i1= i (lower level);
# i2= k (upper level); i3= Z; i5= ionN

# XSTAR data type: 56

const ElectronImpact1Desc = "tabulated collision strength, bautista"

struct ElectronImpact1{I} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    N::I
    T::Tuple
    Υ::Tuple
end

function ElectronImpact1(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ElectronImpact1(Int8(rate), label, ivec..., Tuple(rvec[1:end÷2]),
        Tuple(rvec[end÷2+1:end]))
end

function rate(coef::ElectronImpact1)

end

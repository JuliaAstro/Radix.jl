# Electron-impact effective collision strengths for the k−i transition of the
# N-electron ion (CHIANTI fit citep{1992:Burgess,1997:Dere}): r1 = ∆E (Ryd);
# r2 = C; r3−rjmax+2 = (Υred( j), j= 1,max) (reduced effective collision
# strength); i1= it (transition type); i2= i; i3= k; i4= ionN

# XSTAR data type: 98

const ElectronImpact2Desc = "chianti2016 collisional rates"

struct ElectronImpact2{I,R} <: AbstractRate
    rate::Int8
    label::String
    t::I
    i::I
    k::I
    N::I
    ΔE::R
    C::R
    Υ::Tuple
end

function ElectronImpact2(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ElectronImpact2(Int8(rate), label, ivec..., rvec[1:2]...,
        Tuple(rvec[3:end]))
end

function rate(coef::ElectronImpact2)

end

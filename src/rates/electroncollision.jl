# Electron-impact effective collision strength for the k−i transition of
# N-electron ion (CHIANTI fit [12, 17]): r1 = ∆E (Ryd); r2 = C; r3−r7 =
# (Υ red( j), j= 1,5) (reduced effective collision strength); i1= it
# (transition type); i2= i (lower level); i3= k (upper level); i4= Z; i5= ionN

# XSTAR data type: 51

const ElectronCollisionDesc = "op and chianti line coll rates"

struct ElectronCollision{I, R} <: AbstractRate
    rate::Int8
    label::String
    type::I
    i::I
    k::I
    Z::I
    N::I
    ΔE::R
    C::R
    Υ::Tuple
end

function ElectronCollision(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ElectronCollision(Int8(rate), label, ivec..., rvec[1:2]..., Tuple(rvec[3:7]))
end

function rate(coef::ElectronCollision)

end

# Collision transition rates from superlevels to spectroscopic levels:
# r1−rjnd = (ne( j), j= 1, jnd); rjnd+1−rjnd+nt = (Te( j), j= 1, jnt);
# rjnd+nt+1−rjnd+nt+nt*nd = ((C( j, j′), j′ = 1, j′ nd), j= 1, jnt) (s−1);
# rjnd+nt+nt*nd+1 = λ(Å); i1= nd; i2= nt; i3= i (lower level);
# i4= k (upper level); i5= Z; i6= ionN

# XSTAR data type: 77

const CollisionSuperDesc = "coll rates from 71"

struct CollisionSuper{I, R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    k::I
    Z::I
    N::I
    λ::R
    ne::Tuple
    Te::Tuple
    A::Tuple
end

function CollisionSuper(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Nd, Nt = ivec[1:2]
    CollisionSuper(Int8(rate), label, ivec[3:6]..., rvec[Nd+Nt+Nd*Nt+1],
        Tuple(rvec[1:Nd]), Tuple(rvec[Nd+1:Nd+Nt]),
        Tuple(rvec[Nd+Nt+1:Nd+Nt+Nd*Nt]))
end

function rate(coef::CollisionSuper)
    
end

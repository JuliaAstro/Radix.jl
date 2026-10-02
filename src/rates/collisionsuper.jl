# Collision transition rates from superlevels to spectroscopic levels:
# r1−rjnd = (ne( j), j= 1, jnd); rjnd+1−rjnd+nt = (Te( j), j= 1, jnt);
# rjnd+nt+1−rjnd+nt+nt*nd = ((C( j, j′), j′ = 1, j′ nd), j= 1, jnt) (s−1);
# rjnd+nt+nt*nd+1 = λ(Å); i1= nd; i2= nt; i3= i (lower level);
# i4= k (upper level); i5= Z; i6= ionN

# XSTAR data type: 77

const CollisionSuperDesc = "coll rates from 71"

struct CollisionSuper{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    λ::R                # wavelength (Å)
    ne_grid::Vector{R}  # electron densities (cm⁻³)
    T_grid::Vector{R}   # temperatures (K)
    C::Matrix{R}        # collision rates (s⁻¹), size (ne, T)
end

function CollisionSuper(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Nd, Nt = ivec[1:2]
    CollisionSuper(Int8(rate), label, Transition(ivec[3], ivec[4]), ivec[5:6]..., rvec[Nd+Nt+Nd*Nt+1],
        Vector(rvec[1:Nd]), Vector(rvec[Nd+1:Nd+Nt]),
        reshape(rvec[Nd+Nt+1:Nd+Nt+Nd*Nt], Nd, Nt))
end

# Collision strengths for Fe XIX [10]: r1 = Υ(k,i); i1= i (lower level);
# i2= k (upper level); i3= Z; i4= ionN

# XSTAR data type: 81

const CollisionFe19Desc = "Bhatia Fe XIX collision strengths"

struct CollisionFe19{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I        # atomic number
    ion::I      # ion index (XSTAR ionN)
    Υ::R  # effective collision strength
end

function CollisionFe19(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionFe19(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., rvec[1])    
end

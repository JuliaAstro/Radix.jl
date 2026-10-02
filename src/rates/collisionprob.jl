# Collisional transition probability Cik for N-electron ion computed by
# quantum defect theory (or hydrogenic): i1= 1; i2= i (lower level);
# i3= k (upper level); i4= Z; i5= ionN

# XSTAR data type: 63

const CollisionProbDesc = "h-like cij, bautista (hlike ion)"

struct CollisionProb{I} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    lower::I  # lower level
    upper::I  # upper level
    Z::I      # atomic number
    ion::I    # ion index (XSTAR ionN)
end

function CollisionProb(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionProb(Int8(rate), label, ivec[2:end]...)
end

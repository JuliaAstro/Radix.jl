# Two-photon radiation rate for (k−i) transition of N-electron ion:
# r1 = A(k,i) (s−1); i1= i (lower level); i2= k (upper level); i3= 1; i4= ionN;
#  s1= transition identifier

# XSTAR data type: 76

const TwoPhotonDecayDesc = "2 photon decay"

struct TwoPhotonDecay{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    ion::I    # ion index (XSTAR ionN)
    A::R      # s⁻¹
end

function TwoPhotonDecay(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    TwoPhotonDecay(Int8(rate), label::String, Transition(ivec[1], ivec[2]), ivec[4], rvec[1])
end

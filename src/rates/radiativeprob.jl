# Radiative transition probability Aki for the k−i transition of N-electron
# ion computed by quantum defect theory (or hydrogenic): r1 = 0.0E + 0;
# i1= i (lower level); i2= k (upper level); i3= Z; i5= ionN

# XSTAR data type: 54

const RadiativeProbDesc = "h-like cij, bautista (hlike ion)"

struct RadiativeProb{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I      # atomic number
    ion::I    # ion index (XSTAR ionN)
    A::R      # s⁻¹
end

function RadiativeProb(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    RadiativeProb(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., rvec[1])
end

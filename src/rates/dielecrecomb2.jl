# 

# XSTAR data type: 8

const DielecRecomb2Desc = "dielectronic recombination: arnaud and raymond"

struct DielecRecomb2{R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    C::Vector{R}  # coefficients
    E::Vector{R}  # energies
end

function DielecRecomb2(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    DielecRecomb2(Int8(rate), label, rvec[1:4], rvec[5:8])
end

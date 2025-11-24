# 

# XSTAR data type: 8

const DielecRecomb2Desc = "dielectronic recombination: arnaud and raymond"

struct DielecRecomb2{R} <: AbstractRate
    rate::Int8
    label::String
    cfe::Tuple{R}
    efe::Tuple{R}
end

function DielecRecomb2(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    DielecRecomb2(Int8(rate), label, rvec[1:4], rvec[5:8])
end

function rate(coef::DielecRecomb2)

end

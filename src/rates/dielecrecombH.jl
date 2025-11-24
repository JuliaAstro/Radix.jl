# Dielectronic recombination rate coefficient of the N-electron recombined
# ion [42]: r1 = a; r2 = b; r3= c; r4 = d; r5 = e; r6= f ; i1= ionN

# XSTAR data type: 22

const DielecRecombHDesc = "dielectronic recombination: storey"

struct DielecRecombH <: AbstractRate
    rate::Int8
    label::String
end

function DielecRecombH(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    DielecRecombH(Int8(rate), label)
end

function rate(coef::DielecRecombH)
end

# Dielectronic recombination rate coefficient of N-electron recombined ion
# [1, 2]: r1 = Adi (cm3 s−1 K3/2); r2 = Bdi; r3= T0 (K); r4 = T1 (K); i1= ionN

# XSTAR data type: 7

const DielecRecomb1Desc = "dielectronic recombination: aldrovandi and pequi"

struct DielecRecomb1{R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    A::R   # cm³ s⁻¹ K³ᐟ²
    B::R
    T0::R  # K
    T1::R  # K
end

function DielecRecomb1(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    DielecRecomb1(Int8(rate), label, rvec[1], rvec[2], rvec[3], rvec[4])
end

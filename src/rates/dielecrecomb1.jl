# Dielectronic recombination rate coefficient of N-electron recombined ion
# [1, 2]: r1 = Adi (cm3 s−1 K3/2); r2 = Bdi; r3= T0 (10⁴ K); r4 = T1 (10⁴ K); i1= ionN

# XSTAR data type: 7

const DielecRecomb1Desc = "dielectronic recombination: aldrovandi and pequi"

struct DielecRecomb1{R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    A::R   # cm³ s⁻¹ K³ᐟ²
    B::R
    T0::R  # 10⁴ K, as ucalc uses it with a temperature in 10⁴ K
    T1::R  # 10⁴ K
end

function DielecRecomb1(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    DielecRecomb1(Int8(rate), label, rvec[1], rvec[2], rvec[3], rvec[4])
end

"""
    rate(coef::DielecRecomb1, cell; index=false)

Dielectronic recombination to the ground level of the next ion from the formula of Aldrovandi and
Péquignot (XSTAR ucalc type 7): `frate = nₑ A e^{-T0/T}(1 + B e^{-T1/T}) T^{-3/2}` with T in 10⁴ K,
`T0` and `T1` in the same units and `A` in cm³ s⁻¹ K^{3/2}.
"""
function rate(coef::DielecRecomb1, cell::Cell; index=false, verbose=false)
    index && return (; init=1, final=0, frate=0., irate=0.)
    T = cell.T
    dirt = Float64(coef.A)*T_unit^-1.5*expo(-Float64(coef.T0)/T)*
        (1 + Float64(coef.B)*expo(-Float64(coef.T1)/T))/(T*sqrt(T))
    (; init=1, final=0, frate=cell.nₑ*dirt, irate=0.)
end

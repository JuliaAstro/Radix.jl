# Total dielectronic recombination rate coefficient of N-electron recombined
# ion [http://amdpp.phys.strath.ac.uk/tamoc/DATA/DR/]: r1−rjmax = (C( j), j= 1,
# jmax)(cm3 s−1 K3/2); rjmax+1−rj2*max = (T( j), j= 1, jmax) (K); i1= Z;
# i2= N−1; i3= M; i4= W, i5= ionN

# XSTAR data type: 39

const TotDielecRecombDesc = "total dr  from badnell amdpp.phys.strath.ac.uk"

struct TotDielecRecomb{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    Z::I                # atomic number
    N::I                # number of electrons minus one
    M::I
    W::I
    ion::I              # ion index (XSTAR ionN)
    C::Vector{R}        # cm³ s⁻¹ K^{3/2}
    T::Vector{R}        # K
end

function TotDielecRecomb(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    n = length(rvec) ÷ 2
    TotDielecRecomb(Int8(rate), label, ivec..., Vector(rvec[1:n]), Vector(rvec[n + 1:2n]))
end

"""
    rate(coef::TotDielecRecomb, cell; index=false)

Total dielectronic recombination to the ground level of the next ion (XSTAR ucalc type 39), the sum
`frate = nₑ T^{-3/2} Σᵢ Cᵢ e^{-Tᵢ/T}` of the terms of Badnell's fit, with T in 10⁴ K and the Tᵢ in K.
"""
function rate(coef::TotDielecRecomb, cell::Cell; index=false)
    index && return (; init=1, final=0, frate=0., irate=0.)
    T = cell.T
    dirt = sum(Float64(c)*exp(-Float64(t)/T_unit/T) for (c, t) in zip(coef.C, coef.T); init=0.0)
    (; init=1, final=0, frate=cell.nₑ*dirt*T_unit^-1.5*T^-1.5, irate=0.)
end

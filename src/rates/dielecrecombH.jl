# Dielectronic recombination rate coefficient of the N-electron recombined
# ion [42]: r1 = a; r2 = b; r3= c; r4 = d; r5 = e; i1= ionN

# XSTAR data type: 22

const DielecRecombHDesc = "dielectronic recombination: storey"

struct DielecRecombH{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    ion::I               # ion index (XSTAR ionN)
    coeffs::Vector{R}    # a, b, c, d, e of the fit
end

function DielecRecombH(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    DielecRecombH(Int8(rate), label, ivec[1], Vector(rvec))
end

const storey_T_max = 6.0          # no rate above this temperature (10⁴ K)
const storey_unit = 1e-12

"""
    rate(coef::DielecRecombH, cell; index=false)

Dielectronic recombination to the ground level of the next ion from the fit of Storey (XSTAR ucalc type
22): `frate = nₑ 10⁻¹² (a/T + b + T(c + T d)) T^{-3/2} e^{-e/T}` with T in 10⁴ K, floored at 0 and zero
above 6×10⁴ K.
"""
function rate(coef::DielecRecombH, cell::Cell; index=false, verbose=false)
    none = (; init=1, final=0, frate=0., irate=0.)
    (index || cell.T > storey_T_max) && return none
    a, b, c, d, e = Float64.(coef.coeffs[1:5])
    T = cell.T
    dirt = storey_unit*(a/T + b + T*(c + T*d))*T^-1.5*expo(-e/T)
    (; init=1, final=0, frate=cell.nₑ*max(dirt, 0.0), irate=0.)
end

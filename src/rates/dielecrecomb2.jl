# Dielectronic recombination rate coefficient of N-electron recombined ion [Arnaud and Raymond]:
# r1−r4 = c(1..4) (cm3 s−1 K3/2); r5−r8 = E(1..4) (eV)

# XSTAR data type: 8

const DielecRecomb2Desc = "dielectronic recombination: arnaud and raymond"

struct DielecRecomb2{R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    C::Vector{R}  # coefficients (cm³ s⁻¹ K^{3/2})
    E::Vector{R}  # energies (eV)
end

function DielecRecomb2(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    DielecRecomb2(Int8(rate), label, rvec[1:4], rvec[5:8])
end

"""
    rate(coef::DielecRecomb2, cell; index=false)

Dielectronic recombination to the ground level of the next ion from the sum of four terms of Arnaud and Raymond
(XSTAR ucalc type 8): `frate = nₑ T^{-3/2} Σᵢ Cᵢ e^{-Eᵢ/kT}` with T in 10⁴ K and the energies in eV. The
database has no records of this type: it is checked only against `ucalc` with synthetic records
(`test/reference/ucalc/synthetic.jl`).
"""
function rate(coef::DielecRecomb2, cell::Cell; index=false)
    index && return (; init=1, final=0, frate=0., irate=0.)
    K = constants()
    T = cell.T
    dirt = sum(Float64(c)*expo(-Float64(e)/(K.kT_eV*T)) for (c, e) in zip(coef.C, coef.E); init=0.0)
    (; init=1, final=0, frate=cell.nₑ*dirt*T_unit^-1.5*T^-1.5, irate=0.)
end

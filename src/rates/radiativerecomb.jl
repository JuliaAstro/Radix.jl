# Radiative recombination rate coefficient of N-electron recombined ion [1, 2]:
# r1 = Arad (cm3 s−1); r2 = η; i1= ion_N

# XSTAR data type: 1

const RadRecomDesc = "radiative recombination:  aldrovandi and pequign"

struct RadRecomb{R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    A::R  # cm³ s⁻¹
    η::R
end

function RadRecomb(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    RadRecomb(Int8(rate), label, rvec[1], rvec[2])
end

function rate(coef::RadRecomb, cell::Cell; index=false)
    (; init=1, final=0, frate=(index ? 0. : cell.nₑ*coef.A/cell.T^coef.η), irate=0.)
end

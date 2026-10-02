#

# XSTAR data type: 4

const AtomicLineDesc = "line data radiative: mendosa; raymond and smith"

struct AtomicLine{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    lower::I  # lower level
    upper::I  # upper level
    λ::R      # wavelength
    f::R      # oscillator strength
    A::R      # Einstein A (s⁻¹)
end

function AtomicLine(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AtomicLine(Int8(rate), label, ivec[1], ivec[2], rvec[1], rvec[2], rvec[5])
end

function rate(coef::AtomicLine, cell1::Cell, cell2::Cell)
    if index
        irate, frate = 0., 0.
    else
        init, final = coef.n < coef.m ? (coef.upper, coef.lower) : (coef.lower, coef.upper)
        elin, flin = abs(coef.λ), coef.f
        eeup, eelo = cell1.e, cell2.e
        ggup, gglo = cell1.g, cell2.g
        a = coef.A
        T1, T2 = cell1.T, cell2.T
    end
    (; init=coef.lower, final=coef.upper, )
end

# 

# XSTAR data type: 9

const ChargeExHeDesc = "charge exch. H0 Kingdon and Ferland"

struct ChargeExHe{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    level::I         # level index
    parent_level::I  # level of the recombined ion
    a::R
    b::R
    c::R
    d::R
    T1::R            # K
    T2::R            # K
    ΔE::R            # ΔE/k (10⁴ K)
end

function ChargeExHe(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    # a single index means "ground level" (XSTAR: idest1=1, idest2=nlevp)
    iv = length(ivec) > 1 ? ivec[1:2] : Int32[1, 1]
    ChargeExHe(Int8(rate), label, iv..., rvec...)
end

function rate(coef::ChargeExHe, cell::Cell; index=false, verbose=false,
    ndit=1, nlev=1)
    T = cell.T
    res = 1e-9*coef.a*min(T, 1000.0)^coef.b*(1 + coef.c*expo(coef.d*T))

    init, final, frate = ndit > 1 ? (coef.level, nlev + coef.parent_level - 1, res/6) :
        (1, nlev, 0.0)

    (; init=init, final=final, frate=frate, irate=index ? 0. : res)
end

# Total radiative recombination rate (hydrogenic) for N-electron recombined ion
# [21]: i1= Z; i2= ionN

# XSTAR data type: 30

const TotRadRecombHDesc = "radiative recombination hydrogenic: gould and t"

const C30a = 6.34/100
const C30b = 0.06376
const C30c = 0.9
const C30d = 1.735
const C30e = -1.202
const C30f = -0.298
const C30g = 0.2525
const C30h = 2.105e-22
const C30i = 3.10782e7

struct TotRadRecombH{I} <: AbstractRate
    rate::Int8
    label::String
    Z::I
    N::I
end

function TotRadRecombH(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    TotRadRecombH(Int8(rate), label, ivec...)
end

function rate(coef::TotRadRecombH, cell::Cell; index=false, verbose=false)
    if index
        res = (; init=1, final=0, frate=0., rrate=0.)
    else
        β = coef.nmax^2/(C30a*cell.T)
        # fudge factor makes the two expressions join smoothly
        ypow = min(1, C30b/β^2)
        ϕ1 = (C30c*(1-ypow) + 2/3*ypow)*(C30d + log(β) + 1/6/β)/2
        ϕ2 = β*(C30e*log(β) + C30f)
        ϕ = β < C30g : ϕ2 : ϕ1
        frate = cell.nₑ*2*C30h*C30i*sqrt(cell.T)*β*ϕ
        res = (init=1, final=0, frate=frate, rrate=0.)
    end
    res
end

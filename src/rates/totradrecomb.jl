# Total radiative recombination rate coefficient of N-electron recombined ion
# [http://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/]: r1 = A(cm3 s−1); r2 = B;
# r3 = T0 (K); r4 = T1 (K); r5 = C; r6= T2 (K); i1= Z; i2= N−1; i3= M; i4= W;
# i5= ionN

# XSTAR data type: 38

const TotRadRecomDesc = "total rr  from badnell amdpp.phys.strath.ac.uk"

struct TotRadRecomb{I, R} <: AbstractRate
    rate::Int8
    label::String
    Z::I
    Nm::I
    M::I
    W::I
    N::I
    A::R
    B::R
    T0::R
    T1::R
    C::R
    T2::R
end

function TotRadRecomb(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}
    if length(rvec) > 4
        TotRadRecomb(Int8(rate), label, ivec..., rvec...)
    else
        TotRadRecomb(Int8(rate), label, ivec..., rvec..., 0f0, 0f0)
    end
end

function rate(coef::TotRadRecomb, cell::Cell; index=false, verbose=false)

    if index
        res = (; init=1, final=0, rate=0.)
    else
        b = coef.b + coef.c*exp(-coef.T2/cell.T/1e4)
        term1 = sqrt(cell.T/coef.T0)
        term2 = (1.0 + sqrt(cell.T/coef.T0))^(1.0 - b)
        term3 = (1.0 + sqrt(cell.T/coef.T1))^(1.0 + b)
        frate = cell.nₑ*coef.a/(1e-48 + term1*term2*term3)
        if verbose println() end
        res = (; init=1, final=0, frate=frate)
    end
    res
end

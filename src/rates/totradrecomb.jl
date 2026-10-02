# Total radiative recombination rate coefficient of N-electron recombined ion
# [http://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/]: r1 = A(cm3 s−1); r2 = B;
# r3 = T0 (K); r4 = T1 (K); r5 = C; r6= T2 (K); i1= Z; i2= N−1; i3= M; i4= W;
# i5= ionN

# XSTAR data type: 38

const TotRadRecomDesc = "total rr  from badnell amdpp.phys.strath.ac.uk"

struct TotRadRecomb{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    Z::I                 # atomic number
    parent_electrons::I  # electrons of the recombining ion (N-1)
    M::I
    W::I
    ion::I               # ion index (XSTAR ionN)
    A::R                 # cm³ s⁻¹
    B::R
    T0::R                # K
    T1::R                # K
    C::R
    T2::R                # K
end

function TotRadRecomb(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}
    if length(rvec) > 4
        TotRadRecomb(Int8(rate), label, ivec..., rvec...)
    else
        TotRadRecomb(Int8(rate), label, ivec..., rvec..., 0f0, 0f0)
    end
end

"""
    rate(coef::TotRadRecomb, cell; index=false)

Total radiative recombination (XSTAR ucalc type 38, Badnell's fit). The
characteristic temperatures `T0`, `T1`, `T2` are in K while `cell.T` is in 10⁴ K.
"""
function rate(coef::TotRadRecomb, cell::Cell; index=false, verbose=false)
    index && return (; init=1, final=0, frate=0., irate=0.)
    T = cell.T
    T0, T1 = coef.T0/1e4, coef.T1/1e4
    b = coef.B + coef.C*exp(-coef.T2/1e4/T)
    term1 = sqrt(T/T0)
    term2 = (1.0 + sqrt(T/T0))^(1.0 - b)
    term3 = (1.0 + sqrt(T/T1))^(1.0 + b)
    (; init=1, final=0, frate=cell.nₑ*coef.A/(1e-48 + term1*term2*term3), irate=0.)
end

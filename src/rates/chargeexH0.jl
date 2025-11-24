# H0 charge exchange rate coefficient of N-electron recombined ion [34]:
# r1 = a (10−9 cm3 s−1); r2 = b; r3= c; r4 = d; r5 = T1 (K); r6= T2 (K);
# r7 = ∆E/k (104 K); i1= ionN ; s1= recombining ion identifier

# XSTAR data type: 2

const ChargeExH0Desc = "charge exch. h0: Kingdon and Ferland"

struct ChargeExH0{I, R} <: AbstractRate
    rate::Int8
    label::String
    i::I
    a::R
    b::R
    c::R
    d::R
    T1::R
    T2::R
    ΔE::R

end

function ChargeExH0(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ChargeExH0(Int8(rate), label, ivec[1], rvec...)
end

function rate(coef::ChargeExH0, n::Integer, cell::Cell; index=false, verbose=false)

    if index || T > 5
        res = (; init=1, final=n, frate=0., irate=0.)
    else
        rate = 1e-9*cell.nₕ*coef.a*expo(coef.b*log(cell.T)) * 
            max(0, (1 + coef.c*expo(coef.d*cell.T)))
        frate, irate = coef.rate == 5 ? (0., rate) : (rate, 0.)
        res = (; init=1, final=n, frate=frate, irate=irate)
    end
    if verbose println() end
    res
end

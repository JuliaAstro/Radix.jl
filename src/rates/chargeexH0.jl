# H0 charge exchange rate coefficient of N-electron recombined ion [34]:
# r1 = a (10−9 cm3 s−1); r2 = b; r3= c; r4 = d; r5 = T1 (K); r6= T2 (K);
# r7 = ∆E/k (104 K); i1= ionN ; s1= recombining ion identifier

# XSTAR data type: 2

const ChargeExH0Desc = "charge exch. h0: Kingdon and Ferland"
const cxH0_Tmax = 5.0                 # no rate above this temperature (10⁴ K)

struct ChargeExH0{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    ion::I  # ion index (XSTAR ionN)
    a::R
    b::R
    c::R
    d::R
    T1::R   # K
    T2::R   # K
    ΔE::R   # ΔE/k (10⁴ K)
end

function ChargeExH0(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ChargeExH0(Int8(rate), label, ivec[1], rvec...)
end

function rate(coef::ChargeExH0, cell::Cell; index=false, nlev=0)

    if index || cell.T > cxH0_Tmax
        res = (; init=1, final=nlev, frate=0., irate=0.)
    else
        rate = cx_unit*cell.nₕ*coef.a*expo(coef.b*log(cell.T)) * 
            max(0, (1 + coef.c*expo(coef.d*cell.T)))
        frate, irate = coef.rtype == 5 ? (0., rate) : (rate, 0.)
        res = (; init=1, final=nlev, frate=frate, irate=irate)
    end
    res
end

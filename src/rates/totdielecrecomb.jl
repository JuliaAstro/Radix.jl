# Total dielectronic recombination rate coefficient of N-electron recombined
# ion [http://amdpp.phys.strath.ac.uk/tamoc/DATA/DR/]: r1−rjmax = (C( j), j= 1,
# jmax)(cm3 s−1 K3/2); rjmax+1−rj2*max = (T( j), j= 1, jmax) (K); i1= Z;
# i2= N−1; i3= M; i4= W, i5= ionN

# XSTAR data type: 39

const TotDielecRecombDesc = "total dr  from badnell amdpp.phys.strath.ac.uk"

struct TotDielecRecomb <: AbstractRate
    rate::Int8
    label::String
end

function TotDielecRecomb(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    TotDielecRecomb(Int8(rate), label)
end

function rate(coef::TotDielecRecomb)
end

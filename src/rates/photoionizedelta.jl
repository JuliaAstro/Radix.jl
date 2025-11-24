# Delta functions to add to photoionization cross sections to match ADF DR
# rates: r1 = E(∞) (eV); r1−rjm = (E( j), j= 1, jm) (eV); rjm+1−rj2m = 
# ( f ( j), j= 1, jm) (cm2); i1= n; i2= L; i3= 2S + 1; i4= Z; i5= kN−1;
# i6= ionN−1; i7= iN ; i8= ionN

# XSTAR data type: 74

const PhotoionizeDeltaDesc = "Delta functions to add to phot. x-sections  DR"

struct PhotoionizeDelta{I,R} <: AbstractRate
    rate::Int8
    label::String
    n::I
    L::I
    S2::I
    Z::I
    k::I
    Nm::I
    i::I
    N::I
    E∞::R
    E::Tuple
    f::Tuple
end

function PhotoionizeDelta(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    N = (length(rvec)-1)÷2
    PhotoionizeDelta(Int8(rate), label, ivec..., rvec[1], Tuple(rvec[2:N]),
        Tuple(rvec[N+1:2*N]))
end

function rate(coef::PhotoionizeDelta)
    
end

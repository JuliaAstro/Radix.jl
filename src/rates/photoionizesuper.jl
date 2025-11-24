# Coefficients for recombination and photoionization cross sections of
# superlevels: r1−rjnd = (ne( j), j= 1, jnd); rjnd+1−rjnd+nt = (Te( j),
# j= 1, jnt); rjnd+nt+1−rjnd+nt+nt*nd = ((log α( j, j′), j′ = 1, j′ nd),
# j= 1, jnt); rjnd+nt+nt*nd+1−rjnd+nt+nt*nd+2*nx = (E( j),σ( j), j= 1, jnx);
# i1= nd; i2= nt; i3= nx; i4= n; i5= L; i6= 2S + 1; i7= Z; i8= kN−1;
# i9= ionN−1; i10= iN ; i11= ionN

# XSTAR data type: 70

const PhotoionizeSuperDesc = "Coefficients for phot x-section of suplevels"

struct PhotoionizeSuper{I} <: AbstractRate
    rate::Int8
    label::String
    n::I
    L::I
    S2::I
    Z::I
    kNm::I
    Nm::I
    i::I
    N::I
    ne::Tuple
    Te::Tuple
    α::Tuple
    E::Tuple
    σ::Tuple
end

function PhotoionizeSuper(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Nd, Nt, Nx = ivec[1:3]
    PhotoionizeSuper(Int8(rate), label, ivec[4:end]..., Tuple(rvec[1:Nd]),
        Tuple(rvec[Nd+1:Nd+Nt]), Tuple(rvec[Nd+Nt+1:Nd+Nt+Nd*Nt]),
        Tuple(rvec[Nd+Nt+Nd*Nt+1:2:end-1]), Tuple(rvec[Nd+Nt+Nd*Nt+2:2:end]))
end

function rate(coef::PhotoionizeSuper)

end

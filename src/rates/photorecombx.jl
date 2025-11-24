# Coefficients for recombination and photoionization cross sections of
# superlevels: r1−rjnd = (ne( j), j=1, jnd); rjnd+1−rjnd+nt = (Te( j), j= 1,
# jnt); rjnd+nt+1−rjnd+nt+nt*nd = ((α( j, j′), j′= 1, j′nd), j= 1, jnt);
# rjnd+nt+nt*nd+1−rjnd+nt+nt*nd+2*nx = (E( j),σ( j), j= 1, jnx); i1= nd;
# i2= nt; i3= nx; i4= n; i5= L; i6= 2S + 1; i7= Z; i8= kN−1; i9= ionN−1;
# i10= iN ; i11= ionN

# XSTAR data type: 99

const PhotoRecombXDesc = ""

struct PhotoRecombX{I} <: AbstractRate
    rate::Int8
    label::String
    nd::I
    nt::I
    nx::I
    n::I
    L::I
    S2p::I
    Z::I
    kNm::I
    ionNm::I
    iN::I
    ionN::I
    ne::Tuple                 # electron density
    Te::Tuple                 # electron temperature
    α::Tuple                  # 
    E::Tuple                  # energy
    σ::Tuple                  # cross-section
end

function PhotoRecombX(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    nd, nt, nx = ivec[1:3]
    PhotoRecombX(Int8(rate), label, ivec..., Tuple(rvec[1:nd]),
        Tuple(rvec[nd+1:nd+nt]), Tuple(rvec[nd+nt+1:nd+nt+nd*nt]),
        Tuple(rvec[nd+nt+nd*nt+1:2:nd+nt+nd*nt+2*nx-1]),
        Tuple(rvec[nd+nt+nd*nt+2:2:nd+nt+nd*nt+2*nx]))
end

function rate(coef::PhotoRecombX)
    
end

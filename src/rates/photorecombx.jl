# Coefficients for recombination and photoionization cross sections of
# superlevels: r1−rjnd = (ne( j), j=1, jnd); rjnd+1−rjnd+nt = (Te( j), j= 1,
# jnt); rjnd+nt+1−rjnd+nt+nt*nd = ((α( j, j′), j′= 1, j′nd), j= 1, jnt);
# rjnd+nt+nt*nd+1−rjnd+nt+nt*nd+2*nx = (E( j),σ( j), j= 1, jnx); i1= nd;
# i2= nt; i3= nx; i4= n; i5= L; i6= 2S + 1; i7= Z; i8= kN−1; i9= ionN−1;
# i10= iN ; i11= ionN

# XSTAR data type: 99

const PhotoRecombXDesc = ""

struct PhotoRecombX{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I                # principal quantum number
    L::I                # orbital angular momentum
    spin_mult::I        # 2S+1
    Z::I                # atomic number
    parent::Parent{I}           # parent ion and level
    level::I            # level index
    ion::I              # ion index (XSTAR ionN)
    ne_grid::Vector{R}  # electron densities (cm⁻³)
    T_grid::Vector{R}   # temperatures (K)
    α::Matrix{R}        # α, size (ne, T)
    E_grid::Vector{R}   # energies
    σ::Vector{R}        # cross sections
end

function PhotoRecombX(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    nd, nt, nx = ivec[1:3]
    PhotoRecombX(Int8(rate), label, ivec[4:7]..., Parent(ivec[9], ivec[8]), ivec[10:11]..., Vector(rvec[1:nd]),
        Vector(rvec[nd+1:nd+nt]), reshape(rvec[nd+nt+1:nd+nt+nd*nt], nd, nt),
        Vector(rvec[nd+nt+nd*nt+1:2:nd+nt+nd*nt+2*nx-1]),
        Vector(rvec[nd+nt+nd*nt+2:2:nd+nt+nd*nt+2*nx]))
end

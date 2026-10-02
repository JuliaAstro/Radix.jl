# Coefficients for recombination and photoionization cross sections of
# superlevels: r1−rjnd = (ne( j), j= 1, jnd); rjnd+1−rjnd+nt = (Te( j),
# j= 1, jnt); rjnd+nt+1−rjnd+nt+nt*nd = ((log α( j, j′), j′ = 1, j′ nd),
# j= 1, jnt); rjnd+nt+nt*nd+1−rjnd+nt+nt*nd+2*nx = (E( j),σ( j), j= 1, jnx);
# i1= nd; i2= nt; i3= nx; i4= n; i5= L; i6= 2S + 1; i7= Z; i8= kN−1;
# i9= ionN−1; i10= iN ; i11= ionN

# XSTAR data type: 70

const PhotoionizeSuperDesc = "Coefficients for phot x-section of suplevels"

struct PhotoionizeSuper{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I                # principal quantum number
    L::I                # orbital angular momentum
    spin_mult::I        # 2S+1
    Z::I                # atomic number
    parent_level::I     # level of the parent (N-1 electron) ion
    parent_ion::I       # parent (N-1 electron) ion index
    level::I            # level index
    ion::I              # ion index (XSTAR ionN)
    ne_grid::Vector{R}  # electron densities (cm⁻³)
    T_grid::Vector{R}   # temperatures (K)
    logα::Matrix{R}     # log α, size (ne, T)
    E_grid::Vector{R}   # energies
    σ::Vector{R}        # cross sections
end

function PhotoionizeSuper(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Nd, Nt, Nx = ivec[1:3]
    PhotoionizeSuper(Int8(rate), label, ivec[4:end]..., Vector(rvec[1:Nd]),
        Vector(rvec[Nd+1:Nd+Nt]), reshape(rvec[Nd+Nt+1:Nd+Nt+Nd*Nt], Nd, Nt),
        Vector(rvec[Nd+Nt+Nd*Nt+1:2:end-1]), Vector(rvec[Nd+Nt+Nd*Nt+2:2:end]))
end

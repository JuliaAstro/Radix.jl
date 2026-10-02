# Delta functions to add to photoionization cross sections to match ADF DR
# rates: r1 = E(∞) (eV); r1−rjm = (E( j), j= 1, jm) (eV); rjm+1−rj2m = 
# ( f ( j), j= 1, jm) (cm2); i1= n; i2= L; i3= 2S + 1; i4= Z; i5= kN−1;
# i6= ionN−1; i7= iN ; i8= ionN

# XSTAR data type: 74

const PhotoionizeDeltaDesc = "Delta functions to add to phot. x-sections  DR"

struct PhotoionizeDelta{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I               # principal quantum number
    L::I               # orbital angular momentum
    spin_mult::I       # 2S+1
    Z::I               # atomic number
    parent::Parent{I}           # parent ion and level
    level::I           # level index
    ion::I             # ion index (XSTAR ionN)
    E_inf::R           # eV
    E_grid::Vector{R}  # energies (eV)
    f::Vector{R}       # cm²
end

function PhotoionizeDelta(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    N = (length(rvec)-1)÷2
    PhotoionizeDelta(Int8(rate), label, ivec[1:4]..., Parent(ivec[6], ivec[5]), ivec[7:8]..., rvec[1], Vector(rvec[2:N]),
        Vector(rvec[N+1:2*N]))
end

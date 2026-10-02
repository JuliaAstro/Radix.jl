# Photoionization cross section damped excess of iN -th level of the N-electron
# ion leaving the (N−1)-electron ion in superlevel_[K] kN−1: r1−rjmax = 
# (E( j),σ(E( j)), j= 1, jmax) (Energy in Ryd relative to E(∞), cross section
# in Mb); i1= n; i2= L; i3= 2J; i4= Z; i5= kN−1; i6= iN ; i7= ionN.

# XSTAR data type: 88

const PhotoionizeDampDesc = "Iron inner shell resonance excitation (Patrick)"

struct PhotoionizeDamp{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I               # principal quantum number
    L::I               # orbital angular momentum
    twoJ::I            # 2J
    Z::I               # atomic number
    parent_level::I    # level of the parent (N-1 electron) ion
    superlevel::I      # superlevel index (0 if absent)
    level::I           # level index
    ion::I             # ion index (XSTAR ionN)
    E_grid::Vector{R}  # energies (Ry)
    σ::Vector{R}       # cross sections (Mb)
end

function PhotoionizeDamp(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    iv = length(ivec) < 8 ? (ivec[1:5]..., Int32(0), ivec[6:7]...) : Vector(ivec)
    PhotoionizeDamp(Int8(rate), label, iv..., Vector(rvec[1:2:end-1]),
        Vector(rvec[2:2:end]))
end

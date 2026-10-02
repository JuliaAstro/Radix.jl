# TOPbase partial photoionization cross section (resonance averaged) of iN-th
# level of the N-electron ion leaving the (N−1)-electron ion in the kN−1-th
# level: r1−rj2*max = (E( j),σ(E( j)), j= 1, jmax) (Energy in Ryd relative to
# E(∞), cross section in Mb); i1= n; i2= L; i3= 2J; i4= Z; i5= kN−1;
# i6= ionN−1; i7= iN ; i8= ionN

# XSTAR data type: 53

const ParPhotoIonize2Desc = "op pi xsections"

struct ParPhotoIonize2{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I               # principal quantum number
    L::I               # orbital angular momentum
    twoJ::I            # 2J
    Z::I               # atomic number
    parent_level::I    # level of the parent (N-1 electron) ion
    parent_ion::I      # parent (N-1 electron) ion index
    level::I           # level index
    ion::I             # ion index (XSTAR ionN)
    E_grid::Vector{R}  # energies (Ry)
    σ::Vector{R}       # cross sections (Mb)
end

function ParPhotoIonize2(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ParPhotoIonize2(Int8(rate), label, ivec..., Vector(rvec[1:2:end-1]),
        Vector(rvec[2:2:end]))
end

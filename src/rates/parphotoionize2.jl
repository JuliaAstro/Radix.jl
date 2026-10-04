# TOPbase partial photoionization cross section (resonance averaged) of iN-th
# level of the N-electron ion leaving the (N−1)-electron ion in the kN−1-th
# level: r1−rj2*max = (E( j),σ(E( j)), j= 1, jmax) (Energy in Ryd relative to
# E(∞), cross section in Mb); i1= n; i2= L; i3= 2J; i4= Z; i5= kN−1;
# i6= ionN−1; i7= iN ; i8= ionN

# XSTAR data type: 53

const ParPhotoIonize2Desc = "op pi xsections"

struct ParPhotoIonize2{I, R, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I               # principal quantum number
    L::I               # orbital angular momentum
    twoJ::I            # 2J
    Z::I               # atomic number
    parent::Parent{I}           # parent ion and level
    level::I           # level index
    ion::I             # ion index (XSTAR ionN)
    E_grid::Vector{R}  # energies (Ry)
    σ::Vector{R}       # cross sections (Mb)
    levels::L              # the level data of the database (a Levels)
end

function ParPhotoIonize2(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ParPhotoIonize2(Int8(rate), label, ivec[1:4]..., Parent(ivec[6], ivec[5]), ivec[7:8]..., Vector(rvec[1:2:end-1]),
        Vector(rvec[2:2:end]), levels)
end

"""
    rate(coef::ParPhotoIonize2, cell; levels, radiation, nlev, ptmp=(0.5, 0.5), abund=(0, 0), lfast=1, opacity=nothing)

Photoionization from TOPbase's resonance-averaged cross sections (XSTAR ucalc type
53): as `ParPhotoIonize1` but without the E⁻³ tail, and the threshold includes the
excitation energy of an excited parent level. See `photoionize_level`.
"""
rate(coef::ParPhotoIonize2, cell::Cell; kw...) =
    photoionize_level(coef, cell; extrapolate=false, shifted=true, kw...)

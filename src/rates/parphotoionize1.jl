# Partial photoionization cross section of iN -th level of the N-electron ion
# leaving the (N−1)-electron ion in the kN−1-th level: r1−rj2*max = 
# (E( j),σ(E( j)), j= 1, jmax) (Energy in Ryd relative to E(∞), cross section
# in Mb); i1= n; i2= L; i3= 2J; i4= Z; i5= kN−1; i6= ionN−1; i7= iN ; i8= ionN

# XSTAR data type: 49

const ParPhotoIonize1Desc = "op pi xsections for inner shells"

struct ParPhotoIonize1{I, R, L<:LevelTable} <: AbstractRate
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
    levels::L              # the level data of the database (a LevelTable)
end

function ParPhotoIonize1(rate::Int32, label::String, ivec::I, rvec::R, levels::LevelTable) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ParPhotoIonize1(Int8(rate), label, ivec[1:4]..., Parent(ivec[6], ivec[5]), ivec[7:8]..., Vector(rvec[1:2:end-1]),
        Vector(rvec[2:2:end]), levels)
end

"""
    rate(coef::ParPhotoIonize1, cell; levels, radiation, nlev, ptmp=(0.5, 0.5), abund=(0, 0), lfast=1, opacity=nothing)

Photoionization of a level (XSTAR ucalc type 49) from a cross-section table that is
extended with an E⁻³ tail; see `photoionize_level` for the arguments and the result.
"""
rate(coef::ParPhotoIonize1, cell::Cell; kw...) =
    photoionize_level(coef, cell; extrapolate=true, shifted=false, kw...)

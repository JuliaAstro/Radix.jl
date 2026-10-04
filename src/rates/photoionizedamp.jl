# Photoionization cross section damped excess of iN -th level of the N-electron
# ion leaving the (N−1)-electron ion in superlevel_[K] kN−1: r1−rjmax = 
# (E( j),σ(E( j)), j= 1, jmax) (Energy in Ryd relative to E(∞), cross section
# in Mb); i1= n; i2= L; i3= 2J; i4= Z; i5= kN−1; i6= iN ; i7= ionN.

# XSTAR data type: 88

const PhotoionizeDampDesc = "Iron inner shell resonance excitation (Patrick)"

@with_levels struct PhotoionizeDamp{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I               # principal quantum number
    L::I               # orbital angular momentum
    twoJ::I            # 2J
    Z::I               # atomic number
    parent::Parent{I}           # parent ion (0: not stored) and level
    superlevel::I      # superlevel index (0 if absent)
    level::I           # level index
    ion::I             # ion index (XSTAR ionN)
    E_grid::Vector{R}  # energies (Ry)
    σ::Vector{R}       # cross sections (Mb)
end

function PhotoionizeDamp(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    iv = length(ivec) < 8 ? (ivec[1:5]..., Int32(0), ivec[6:7]...) : Vector(ivec)
    PhotoionizeDamp(Int8(rate), label, iv[1:4]...,
        Parent(zero(eltype(iv)), iv[5]), iv[6:8]..., Vector(rvec[1:2:end-1]),
        Vector(rvec[2:2:end]))
end

"""
    rate(coef::PhotoionizeDamp, cell; levels, radiation, nlev, ptmp=(0.5, 0.5), abund=(0, 0), lfast=1, opacity=nothing)

Photoionization with the damped-excess cross sections (XSTAR ucalc type 88): as
`ParPhotoIonize1` but always to the continuum, and only the photoionization rate
`frate` and the opacity are returned (the recombination and energy terms are zero, as
in ucalc). See `photoionize_level`.
"""
rate(coef::PhotoionizeDamp, cell::Cell; kw...) =
    photoionize_level(coef, cell; extrapolate=true, shifted=false, parent=false,
        rates_only=true, kw...)

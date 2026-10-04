# Fits to fine-structure collision strengths for He-like ions [32]: r1−rjmax =
# coefficients; i1= i (lower level); i2= k (upper level); i3= Z; i4= ionN

# XSTAR data type: 66

const CollisionHeFineDesc = "Like type 69 but, data in fine structure"

struct CollisionHeFine{I, R, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I               # atomic number
    ion::I             # ion index (XSTAR ionN)
    coeffs::Vector{R}  # fit coefficients
    levels::L              # the level data of the database (a Levels)
end

function CollisionHeFine(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionHeFine(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., Vector(rvec), levels)
end

"""
    rate(coef::CollisionHeFine, cell; index=false)

Collisional excitation and de-excitation of a He-like ion from the fine-structure fits (XSTAR ucalc type 66,
`calt66`): the first coefficient is the energy (eV) of the transition, used for the Boltzmann factor, the
temperature floor ΔE/50k and the energy of the rates instead of the energy difference of the levels, which only
order them. Υ is the first term of type 69's fit (and two more if there are more than 6 coefficients, a form
the database does not have) and is not clamped at 0; otherwise like `rate(::CollisionLS, ...)`.
"""
rate(coef::CollisionHeFine, cell::Cell; index=false) =
    helike_fit_collision(coef, cell, index, :fine)

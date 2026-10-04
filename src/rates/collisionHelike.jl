# Analytic fits for effective collision strengths in He-like ions [53]:
# r1−rjmax = coefficients; i1= i (lower level); i2= k (upper level); i3= Z;
# i8= ionN

# XSTAR data type: 68

const CollisionHelikeDesc = ""

struct CollisionHelike{I, R, L<:LevelTable} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I               # atomic number
    ion::I             # ion index (XSTAR ionN)
    coeffs::Vector{R}  # fit coefficients
    levels::L              # the level data of the database (a LevelTable)
end

function CollisionHelike(rate::Int32, label::String, ivec::I, rvec::R, levels::LevelTable) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionHelike(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., Vector(rvec), levels)
end

"""
    rate(coef::CollisionHelike, cell; index=false)

Collisional excitation and de-excitation of a He-like ion from a quadratic fit of the effective collision
strength in log₁₀(T/Z³) (XSTAR ucalc type 68, `calt68`), clamped at 0, with a temperature of at least
ΔE/50k; otherwise like `rate(::CollisionLS, ...)`.
"""
rate(coef::CollisionHelike, cell::Cell; index=false) =
    helike_fit_collision(coef, cell, index, :helike)

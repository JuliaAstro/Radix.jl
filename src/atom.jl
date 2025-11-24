# Atom type

# XSTAR data type: 13

const AtomDesc = "element data:"

"""
Atom(label, i1-i8, mass, Eion)
"""
struct Atom{I, F}
    # rate::Int8 = 13
    label::String
    N::I     # electron number
    Z::I     # atomic number
    mass::F  # atomic mass
    E::F     # atomic ionization energy
end

function Atom(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Atom(label, ivec..., rvec...)
end

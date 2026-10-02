# Atom type

# XSTAR data type: 13

const AtomDesc = "element data:"

"""
Atom(label, i1-i8, mass, Eion)
"""
struct Atom{I, F}
    # rtype = 13
    label::String
    n_ions::I       # number of ionization stages
    Z::I            # atomic number
    abundance::F    # abundance relative to hydrogen
    mass::F         # atomic mass
end

function Atom(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Atom(label, ivec..., rvec...)
end

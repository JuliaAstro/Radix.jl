# Ion type

# XSTAR data type: 14

const IonDesc = "ion data"

"""
Ion(label, Zp, Z, N, E)
"""
struct Ion{I, F}
    # rate::Int8 = 14
    label::String
    I::I            # ionization number
    Z::I            # atomic number
    N::I            # ion number
    E::F            # ionization energy
end

function Ion(rate::Int32, label::String, ivec::I, rvec::F) where
    {I<:AbstractVector{Int32}, F<:AbstractVector{Float32}}

    Ion(label, ivec..., rvec[1])
end

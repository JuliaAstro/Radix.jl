# Ion type

# XSTAR data type: 14

const IonDesc = "ion data"

"""
Ion(label, Zp, Z, N, E)
"""
struct Ion{I, F}
    # rtype = 14
    label::String
    stage::I        # ionization stage (1 = neutral)
    Z::I            # atomic number
    ion::I          # ion index (what rate records call `ion`)
    E_ion::F        # ionization energy (eV)
end

function Ion(rate::Int32, label::String, ivec::I, rvec::F) where
    {I<:AbstractVector{Int32}, F<:AbstractVector{Float32}}

    Ion(label, ivec..., rvec[1])
end

# Data attributes of the i-th level of N-electron ion: r1 = E(i) (eV);
# r2 = (2J + 1); r3= ν(effective quantum number); r4 = E(∞) (eV); i1= n;
# i2= (2S + 1); i3= L; i4= Z; i5= i; i6= ionN ; s1= level configuration
# assignment

# XSTAR data type: 6

const AtomicLevelDesc = "level data"

struct AtomicLevel{I, R} <: AbstractRate
    rate::Int8                # rate type
    label::String             # label
    n::I                      # level
    S2p::I                    # 
    L::I                      # angular momentum
    Z::I                      # atomic number
    i::I                      # initial level
    N::I                      # ionization number
    Ei::R                     # ionization energy
    J2p::R                    # total angular momentum
    ν::R                      # frequency
    Einf::R                   # ionization energy at infinity
end

function AtomicLevel(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AtomicLevel(Int8(rate), label, ivec..., rvec...)
end

function rate(coef::AtomicLevel, cell::Cell;
    index=false, verbose=false)

    (; init=coef.i, final=0, frate=0.)
end

# Data attributes of the i-th level of N-electron ion: r1 = E(i) (eV);
# r2 = (2J + 1); r3= ν(effective quantum number); r4 = E(∞) (eV); i1= n;
# i2= (2S + 1); i3= L; i4= Z; i5= i; i6= ionN ; s1= level configuration
# assignment

# XSTAR data type: 6

const AtomicLevelDesc = "level data"

struct AtomicLevel{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I          # principal quantum number
    spin_mult::I  # 2S+1
    L::I          # orbital angular momentum
    Z::I          # atomic number
    level::I      # level index
    ion::I        # ion index (XSTAR ionN)
    E::R          # level energy (eV)
    g::R          # statistical weight 2J+1
    n_eff::R      # effective quantum number
    E_inf::R      # ionization energy at infinity (eV)
end

function AtomicLevel(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AtomicLevel(Int8(rate), label, ivec..., rvec...)
end

function rate(coef::AtomicLevel, cell::Cell;
    index=false)

    (; init=coef.level, final=0, frate=0., irate=0.)
end

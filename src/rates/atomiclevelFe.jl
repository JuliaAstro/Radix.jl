# Level data for Fe UTA [24]: r1 = E(i) (eV); r2 = (2J +1); r3= 0.0; r4 = 0.0;
# i1= 1; i5= i (level); i6= ionN ; s1= level configuration assignment

# XSTAR data type: 83

const AtomicLevelFeDesc = "Fe UTA level data"

struct AtomicLevelFe{I,R} <: AbstractRate
   rate::Int8                 # rate type
   label::String              # label
   i::I                       # level
   N::I                       # ionization number
   E::R                       # ionization energy
   J2::R                      # 2J value
end

function AtomicLevelFe(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}
    
    AtomicLevelFe(Int8(rate), label, ivec[2], ivec[3], rvec[1], rvec[2])
end

function rate(coef::AtomicLevelFe, cell::Cell)

end

# Collisional ionization rates for N-electron ion [11]: r1 = E(th) (eV);
# r2 = T0 (K); r3−rjmax+2 = (ρ( j), j=1, jmax) (effective collision strength);
# i1= i (level); i5= ionN

# XSTAR data type: 95

const CollisionIonizeDesc = "Bryans CI rates"

struct CollisionIonize{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    level::I      # level index
    ion::I        # ion index (XSTAR ionN)
    E_th::R       # threshold energy (eV)
    T0::R         # K
    ρ::Vector{R}  # effective collision strengths
end

function CollisionIonize(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionIonize(Int8(rate), label, ivec..., rvec[1], rvec[2],
        Vector(rvec[2:end]))
end

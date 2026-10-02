# APED line (k−i) radiation rates [20]: r1 = λ(Å); r2 = 0.0; r3= A(k,i) (s−1);
# i1= i (lower level); i2= k (upper level); i3= Z; i4= ionN

# XSTAR data type: 91

const RadiativeAPEDDesc = "aped line wavelengths same as 50"

struct RadiativeAPED{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    lower::I  # lower level
    upper::I  # upper level
    Z::I      # atomic number
    ion::I    # ion index (XSTAR ionN)
    λ::R      # wavelength (Å)
    A::R      # s⁻¹
end

function RadiativeAPED(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    RadiativeAPED(Int8(rate), label, ivec..., rvec[1], rvec[3])
end

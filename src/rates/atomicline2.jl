# Line (k−i) radiation rates of N-electron ion: r1 = λ(Å); r2 = g f (i,k);
# r3= A(k,i) (s−1); i1= i (lower level); i2= k (upper level); i3= Z; i4= ionN

# XSPEC data type: 50

const AtomicLine2Desc = "op line rad. rates"

struct AtomicLine2{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    lower::I  # lower level
    upper::I  # upper level
    Z::I      # atomic number
    ion::I    # ion index (XSTAR ionN)
    λ::R      # wavelength (Å)
    gf::R     # weighted oscillator strength
    A::R      # Einstein A (s⁻¹)
end

function AtomicLine2(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AtomicLine2(Int8(rate), label::String, ivec..., rvec...)
end

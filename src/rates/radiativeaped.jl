# APED line (k−i) radiation rates [20]: r1 = λ(Å); r2 = 0.0; r3= A(k,i) (s−1);
# i1= i (lower level); i2= k (upper level); i3= Z; i4= ionN

# XSTAR data type: 91

const RadiativeAPEDDesc = "aped line wavelengths same as 50"

struct RadiativeAPED{I, R, L<:LevelTable} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I      # atomic number
    ion::I    # ion index (XSTAR ionN)
    λ::R      # wavelength (Å)
    A::R      # s⁻¹
    levels::L              # the level data of the database (a LevelTable)
end

function RadiativeAPED(rate::Int32, label::String, ivec::I, rvec::R, levels::LevelTable) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    RadiativeAPED(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., rvec[1], rvec[3], levels)
end

"""
    rate(coef::RadiativeAPED, cell; kw...)

Radiative decay and photoexcitation of an APED line (XSTAR ucalc type 91). `ucalc` treats it as type 50
(it jumps to the same code, the wavelength, the stored 0 and A being in the same places), so the
result and the keywords are those of `rate(::AtomicLine2, ...)`.
"""
function rate(coef::RadiativeAPED, cell::Cell; kw...)
    line = AtomicLine2(coef.rtype, "", coef.transition, coef.Z, coef.ion, coef.λ, zero(coef.λ), coef.A, coef.levels)
    rate(line, cell; kw...)
end

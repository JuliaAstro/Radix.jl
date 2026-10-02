# 

# XSTAR data type: 3

const AutoionizeDesc = "AutoIonizeation: hamilton, sarazin chevalier"

const C03a = 0.861707

struct Autoionize{R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    C::R  # rate prefactor
    E::R  # activation energy
end

function Autoionize(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Autoionize(Int8(rate), label, rvec[1], rvec[2])
end

function rate(coef::Autoionize, cell::Cell; index=false, verbose=false)

    frate = index ? 0. : cell.nₑ*coef.C*expo(-coef.E/C03a/cell.T)/sqrt(cell.T)
    (; init=1, final=1, frate=frate, irate=0.)
end

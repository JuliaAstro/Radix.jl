# 

# XSTAR data type: 3

const AutoionizeDesc = "AutoIonizeation: hamilton, sarazin chevalier"

const C03a = 0.861707

struct Autoionize{R} <: AbstractRate
    rate::Int8
    label::String
    cai::R
    eai::R
end

function Autoionize(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Autoionize(Int8(rate), label, rvec[1], rvec[2])
end

function rate(coef::Autoionize, cell::Cell; index=false, verbose=false)

    frate = index ? 0. : den*coef.cai*expo(-coef.eai/(C03a*T))/sqrt(T)
    (; init=1, final=1, frate=frate)
end

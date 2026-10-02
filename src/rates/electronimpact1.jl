# Electron-impact effective collision strengths for the k−i transition of
# N-electron ion: r1−rjmax = (log Te( j), j= 1, jmax) (K); rj(max+1)−rj(2*max)=
# (Υ(Te( j)), j= 1, jmax) (effective collision strength); i1= i (lower level);
# i2= k (upper level); i3= Z; i5= ionN

# XSTAR data type: 56

const ElectronImpact1Desc = "tabulated collision strength, bautista"

struct ElectronImpact1{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # lower and upper level
    Z::I                # atomic number
    ion::I              # ion index (XSTAR ionN)
    T_grid::Vector{R}   # temperatures (K)
    Υ::Vector{R}  # effective collision strengths
end

function ElectronImpact1(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ElectronImpact1(Int8(rate), label, Transition(ivec[1], ivec[2]), ivec[3:4]..., Vector(rvec[1:end÷2]),
        Vector(rvec[end÷2+1:end]))
end

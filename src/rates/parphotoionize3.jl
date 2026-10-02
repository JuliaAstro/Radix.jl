# Partial photoionization cross section of iN -th level of the N-electron ion
# leaving the (N−1)-electron ion in the kN−1-th level [51]: r1 = E(th) (eV);
# r2 = E(0) (eV); r3= σ(0) (Mb); r4 = y(a); r5 = P; r6= y(w) ; i1= N;
# i2= n (shell principal quantum number); i3= l (orbital quantum number of the
# subshell); i4= kN−1; i5= ionN−1; i6= iN ; i7= ionN ; s1= shell-ion identifier

# XSTAR data type: 59

const ParPhotoIonize3Desc = "verner pi x!"

struct ParPhotoIonize3{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n_electrons::I   # electrons in the ion
    n::I             # shell principal quantum number
    l::I             # subshell orbital quantum number
    parent::Parent{I}           # parent ion and level
    level::I         # level index
    ion::I           # ion index (XSTAR ionN)
    E_th::R          # eV
    E0::R            # eV
    σ0::R            # Mb
    ya::R
    P::R
    yw::R
end

function ParPhotoIonize3(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    ParPhotoIonize3(Int8(rate), label, ivec[1:3]..., Parent(ivec[5], ivec[4]), ivec[6:7]..., rvec...)
end

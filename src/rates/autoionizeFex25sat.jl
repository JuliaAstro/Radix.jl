# Autoionization rates for Fe XXIV satellites [8]: r1 = Aa(k,i) (s−1);
# r2 = E(k) (eV above ionization limit); i1= ionN , i2= kN ; i3= ionN−1;
# i4= iN−1; i5= ionN

# XSTAR data type: 75

const AutoionizeFe25SatDesc = "autoionization data for Fe XXiV satellites"

@with_levels struct AutoionizeFe25Sat{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    ion::I           # ion index (XSTAR ionN)
    level::I         # level index
    parent::Parent{I}           # parent ion and level
    A_auto::R        # autoionization rate (s⁻¹)
    E::R             # energy above ionization limit (eV)
end

function AutoionizeFe25Sat(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AutoionizeFe25Sat(Int8(rate), label, ivec[1], ivec[2], Parent(ivec[3], ivec[4]), rvec...)    
end

"""
    rate(coef::AutoionizeFe25Sat, cell; index=false)

Autoionization of an Fe XXIV satellite level (XSTAR ucalc type 75) with the rate of `satellite_capture` (no weight):
`irate` is that rate times the electron density and `frate` is 0. `init` is the level (at least 1) and `final` the
level of the parent ion counted from the continuum (`nlev` + parent level − 1, at least 1).
"""
function rate(coef::AutoionizeFe25Sat, cell::Cell; index=false)
    levels = levels_of(coef)
    nlev = nlevels(levels, coef.ion)
    init = max(Int(coef.level), 1)
    final = max(Int(coef.parent.level) + nlev - 1, 1)
    index && return (; init, final, frate=0., irate=0.)
    irate = satellite_capture(Float64(coef.A_auto), Float64(coef.E), 1.0, cell.T)*cell.nₑ
    (; init, final, frate=0., irate)
end

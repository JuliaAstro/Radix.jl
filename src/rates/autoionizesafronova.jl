# Autoionization rates for satellite levels from Safronova, in Kato et al. (1997, ADNDT 67, 225): r1 = Aa(k,i) (s−1, not used by ucalc);
# r2 = a factor of the rate (the weight); r3 = E(k) (eV above ionization limit); i1= ionN , i2= kN ; i3= ionN−1; i4= iN−1; i5= ionN

# XSTAR data type: 96

const AutoionizeSafronovaDesc = "Autoinization rates (in s^-1) for satellite lvls"

struct AutoionizeSafronova{I, R, L<:Levels} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    ion::I           # ion index (XSTAR ionN)
    level::I         # level index
    parent::Parent{I}           # parent ion and level
    A_auto::R        # autoionization rate (s⁻¹; ucalc does not use it)
    weight::R        # the factor of the rate
    E::R             # energy above ionization limit (eV)
    levels::L              # the level data of the database (a Levels)
end

function AutoionizeSafronova(rate::Int32, label::String, ivec::I, rvec::R, levels::Levels) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AutoionizeSafronova(Int8(rate), label, ivec[1], ivec[2], Parent(ivec[3], ivec[4]), rvec[1:3]..., levels)
end

const safronova_coeff = Float64(2.069f-3)      # the capture rate `2.069e-3 T⁻³ᐟ² e^{-E/kT}` (T in K)

"""
    rate(coef::AutoionizeSafronova, cell; index=false)

Capture into a satellite level (XSTAR ucalc type 96), the inverse of the autoionization of Safronova in Kato et al. (1997): `irate` is
`2.069×10⁻³ T^{-3/2} e^{-E/kT} w nₑ` with the temperature `T` in K, the energy `E` above the ionization limit and the factor `w` (the second real); `frate` is 0 and the
stored autoionization rate is not used. `init` is the level (at least 1) and `final` the level of the parent ion counted from the continuum (`nlev` + parent level − 1, at least 1).
"""
function rate(coef::AutoionizeSafronova, cell::Cell; index=false)
    levels = coef.levels
    nlev = nlevels(levels, coef.ion)
    K = constants()
    init = max(Int(coef.level), 1)
    final = max(Int(coef.parent.level) + nlev - 1, 1)
    index && return (; init, final, frate=0., irate=0.)
    temp = cell.T*T_unit
    irate = safronova_coeff/temp^1.5*expo(-Float64(coef.E)/(K.kT_eV*cell.T))*Float64(coef.weight)*cell.nₑ
    (; init, final, frate=0., irate)
end

# Autoionization rates for Fe XXIV satellites [8]: r1 = Aa(k,i) (s−1);
# r2 = E(k) (eV above ionization limit); i1= ionN , i2= kN ; i3= ionN−1;
# i4= iN−1; i5= ionN

# XSTAR data type: 75

const AutoionizeFe25SatDesc = "autoionization data for Fe XXiV satellites"

struct AutoionizeFe25Sat{I, R} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    ion::I           # ion index (XSTAR ionN)
    level::I         # level index
    parent_ion::I    # parent (N-1 electron) ion index
    parent_level::I  # level of the parent (N-1 electron) ion
    A_auto::R        # autoionization rate (s⁻¹)
    E::R             # energy above ionization limit (eV)
end

function AutoionizeFe25Sat(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    AutoionizeFe25Sat(Int8(rate), label, ivec[1:4]..., rvec...)    
end

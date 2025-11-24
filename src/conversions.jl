# Constants and Conversions

AA2keV(x::Real) = 1e5*hh*ccc/ergsev/x
keV2AA(x::Real) = AA2keV(x)
AA2eV(x::Real)  = 1e3*AA2keV(x)
eV2AA(x::Real)  = keV2AA(1e-3*x)
erg2eV(x::Real) = ergsev*x
eV2erg(x::Real) = x/ergsev

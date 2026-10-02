"""
AbstractRate

All specific rates are subtypes of AbstractRate.

Every rate implements the same interface:

    rate(coef::AbstractRate, cell::Cell; index=false, verbose=false)

and returns a NamedTuple `(; init, final, frate, irate)`, where `init` and
`final` are the levels connected by the transition, `frate` is the forward rate
and `irate` the inverse rate (0 when the rate has no inverse). With
`index=true` only the connectivity is wanted and the rates are returned as 0.
"""
abstract type AbstractRate end

# Fallback for rate types whose physics has not been ported yet.
rate(coef::AbstractRate, cell; kw...) =
    error("rate not implemented for $(typeof(coef))")

# Convenience: coef(cell) is rate(coef, cell), so `coef.(cells)` scans a rate
(coef::AbstractRate)(cell; kw...) = rate(coef, cell; kw...)

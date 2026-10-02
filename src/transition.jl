"""
    Transition(lower, upper)

Lower and upper level indices of a bound-bound transition within one ion.
"""
struct Transition{I<:Integer}
    lower::I
    upper::I
end

"""
    Parent(ion, level)

The parent ion (the N-1 electron ion left behind by photoionization or
autoionization, or recombined from) and the level it is left in. `ion` is 0
when the data do not give the parent ion (`ChargeExHe`, `PhotoionizeDamp`,
`AutoionizeSat`).
"""
struct Parent{I<:Integer}
    ion::I
    level::I
end

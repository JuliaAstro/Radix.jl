# CHIANTI/Burgess & Tully fits to effective collision strengths, as in XSTAR's
# splinem() and upsil().

"""
    splinem(p, x)

5-point spline through the knots `p[1] = y(0)`, `p[2] = y(1/4)`, `p[3] = y(1/2)`,
`p[4] = y(3/4)`, `p[5] = y(1)`, evaluated for `x` in (0, 1).
"""
function splinem(p, x)
    p1, p2, p3, p4, p5 = Float64.(p)
    s  = 1/30
    s2 = 32s*(19p1 - 43p2 + 30p3 - 7p4 + p5)
    s3 = 160s*(-p1 + 7p2 - 12p3 + 7p4 - p5)
    s4 = 32s*(p1 - 7p2 + 30p3 - 43p4 + 19p5)
    if x <= 0.25
        x0 = x - 0.125
        t3 = 0.0
        t2 = 0.5s2
        t1 = 4(p2 - p1)
        t0 = 0.5(p1 + p2) - 0.015625t2
    elseif x <= 0.5
        x0 = x - 0.375
        t3 = 20s*(s3 - s2)
        t2 = 0.25(s2 + s3)
        t1 = 4(p3 - p2) - 0.015625t3
        t0 = 0.5(p2 + p3) - 0.015625t2
    elseif x <= 0.75
        x0 = x - 0.625
        t3 = 20s*(s4 - s3)
        t2 = 0.25(s3 - s4)
        t1 = 4(p4 - p3) - 0.015625t3
        t0 = 0.5(p3 + p4) - 0.015625t2
    else
        x0 = x - 0.875
        t3 = 0.0
        t2 = 0.5s4
        t1 = 4(p5 - p4)
        t0 = 0.5(p4 + p5) - 0.015625t2
    end
    t0 + x0*(t1 + x0*(t2 + x0*t3))
end

"""
    chianti_upsilon(kind, ΔE, C, p, T)

Effective collision strength Υ from the 5-point reduced fit `p` for transition
type `kind` (1–4), transition energy `ΔE` (Ry), abscissa scale `C` and
temperature `T` (K), after Burgess & Tully (1992).
"""
function chianti_upsilon(kind, ΔE, C, p, T)
    e = abs(T/(1.57888e5*ΔE))
    kind in (1, 2, 3, 4) || throw(ArgumentError("unknown transition type $kind"))
    x = kind in (1, 4) ? log((e + C)/C)/log(e + C) : e/(e + C)
    y = splinem(p, x)
    kind == 1 && (y *= log(e + 2.71828))
    kind == 3 && (y /= e + 1)
    kind == 4 && (y *= log(e + C))
    y
end

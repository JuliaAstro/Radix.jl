# CHIANTI/Burgess & Tully fits to effective collision strengths, as in XSTAR's
# splinem(), upsil() and upsiln() (xstarlib/src/*.f90).

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
        t2 = 0.25(s3 + s4)
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

# Natural cubic spline (prepspline/calcspline): second derivatives of the spline
# through (x, y), then evaluation at xt. Outside the nodes the end segments are
# extrapolated, as in XSTAR.
function spline_second_derivatives(x, y)
    n = length(x)
    y2 = zeros(n)
    u = zeros(n)
    for i in 2:n - 1
        sig = (x[i] - x[i-1])/(x[i+1] - x[i-1])
        p = sig*y2[i-1] + 2.0
        y2[i] = (sig - 1.0)/p
        u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1])
        u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p
    end
    y2[n] = 0.0
    for k in n - 1:-1:1
        y2[k] = y2[k]*y2[k+1] + u[k]
    end
    y2
end

function eval_spline(x, y, y2, xt)
    klo, khi = 0, length(x)
    while khi - klo > 1
        k = (khi + klo) ÷ 2
        x[k] > xt ? (khi = k) : (klo = k)
    end
    # below the first node XSTAR reads x(0) (undefined); use the first segment
    klo = max(klo, 1)
    khi = max(khi, klo + 1)
    h = x[khi] - x[klo]
    h == 0 && throw(ArgumentError("calcspline: bad x array"))
    a = (x[khi] - xt)/h
    b = (xt - x[klo])/h
    a*y[klo] + b*y[khi] + ((a^3 - a)*y2[klo] + (b^3 - b)*y2[khi])*h^2/6
end

"""
    chianti_upsilon(kind, ΔE, C, x, Υ, T)

Effective collision strength from an n-point reduced fit (XSTAR `upsiln`) for
transition type `kind` (1–6): `Υ` are the reduced values at the nodes `x` in
[0, 1], `ΔE` the transition energy (Ry), `C` the scale and `T` the temperature (K).
"""
function chianti_upsilon(kind, ΔE, C, x::AbstractVector, Υ::AbstractVector, T)
    kte = T/ΔE/1.57888e5
    kind in 1:6 || throw(ArgumentError("unknown transition type $kind"))
    xt = kind in (1, 4) ? 1 - log(C)/log(kte + C) : kte/(kte + C)
    xs, ys = Float64.(x), Float64.(Υ)
    sups = eval_spline(xs, ys, spline_second_derivatives(xs, ys), xt)
    kind == 1 && return sups*log(kte + exp(1.0))
    kind == 2 && return sups
    kind == 3 && return sups/(kte + 1)
    kind == 4 && return sups*log(kte + C)
    kind == 5 && return sups/kte
    10.0^sups
end

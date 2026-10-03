# CHIANTI/Burgess & Tully fits to effective collision strengths, as in XSTAR's
# splinem(), upsil() and upsiln() (xstarlib/src/*.f90).

# constants of the 5-point spline (splinem) and the fits
const spline_dx = 0.25                  # spacing of the nodes in x
const spline_s = 1/30
const spline_stencils = ((19.0, -43.0, 30.0, -7.0, 1.0),     # second derivative at x = 1/4
                         (-1.0, 7.0, -12.0, 7.0, -1.0),      #   ... at x = 1/2
                         (1.0, -7.0, 30.0, -43.0, 19.0))     #   ... at x = 3/4
const spline_stencil_scale = (32.0, 160.0, 32.0)
const spline_t3_scale = 20.0
const spline_slope = 4.0                # 1/spline_dx
const spline_h2 = 0.015625              # spline_dx²/4
const chianti_Ry_K = 1.57888e5          # K per Ry
const chianti_e = 2.71828               # XSTAR's value of e in the kind-1 factor

"""
    splinem(p, x)

5-point spline through the knots `p[1] = y(0)`, `p[2] = y(1/4)`, `p[3] = y(1/2)`,
`p[4] = y(3/4)`, `p[5] = y(1)`, evaluated for `x` in (0, 1).
"""
function splinem(p, x)
    q = Float64.(p)
    s2, s3, s4 = ntuple(k -> spline_stencil_scale[k]*spline_s*sum(spline_stencils[k] .* q), 3)
    segment = x <= spline_dx ? 1 : x <= 2spline_dx ? 2 : x <= 3spline_dx ? 3 : 4
    x0 = x - (segment - 0.5)*spline_dx
    if segment == 1
        t3 = 0.0
        t2 = s2/2
        t1 = spline_slope*(q[2] - q[1])
        t0 = (q[1] + q[2])/2 - spline_h2*t2
    elseif segment == 2
        t3 = spline_t3_scale*spline_s*(s3 - s2)
        t2 = (s2 + s3)/4
        t1 = spline_slope*(q[3] - q[2]) - spline_h2*t3
        t0 = (q[2] + q[3])/2 - spline_h2*t2
    elseif segment == 3
        t3 = spline_t3_scale*spline_s*(s4 - s3)
        t2 = (s3 + s4)/4
        t1 = spline_slope*(q[4] - q[3]) - spline_h2*t3
        t0 = (q[3] + q[4])/2 - spline_h2*t2
    else
        t3 = 0.0
        t2 = s4/2
        t1 = spline_slope*(q[5] - q[4])
        t0 = (q[4] + q[5])/2 - spline_h2*t2
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
    e = abs(T/(chianti_Ry_K*ΔE))
    kind in (1, 2, 3, 4) || throw(ArgumentError("unknown transition type $kind"))
    x = kind in (1, 4) ? log((e + C)/C)/log(e + C) : e/(e + C)
    y = splinem(p, x)
    kind == 1 && (y *= log(e + chianti_e))
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
    kte = T/ΔE/chianti_Ry_K
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

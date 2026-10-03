# Hydrogenic collision routines used by XSTAR's ucalc type 63 (CollisionProb):
# Gordon's radial integrals (anl1), Pengelly & Seaton / Hummer & Storey l-changing
# rates (velimp), and Johnson / Seaton n-changing excitation (erc, impactn,
# szcoll). Translated from xstarsub.f; the arithmetic is Float64 where the Fortran
# is single precision.

# ---------------------------------------------------------------------------
# Constants of the fits below (local to this file)

# exponential integral (Abramowitz & Stegun 5.1.53 and 5.1.56)
const expint_series = (-0.57721566, 0.99999193, -0.24991055, 0.05519968,
    -0.00976004, 0.00107857)                      # x ≤ 1
const expint_denominator = (9.5733223454, 25.6329561486, 21.0996530827,
    3.9584969228)                                 # x > 1
const expint_numerator = (8.5733287401, 18.0590169730, 8.6347608925,
    0.2677737343)                                 # x > 1

# Gordon (1929) dipole rates

# Pengelly & Seaton (1964) l-changing collisions
const ps_dnl_coeff = 6.0
const ps_rho_coeff = 0.72                         # pa = ps_rho_coeff / total decay rate
const ps_b_coeff = 1.157
const ps_exp_cutoff = 50.0                        # no exponential integral above this
const ps_series_cutoff = 1e-3                     # series instead of closed form below this
const ps_series = (1/3, 1/4, 1/10)                # coefficients of that series

# Seaton (1962) impact-parameter functions
const impcfn_a = (0.9947187, 0.6030883, -2.372843, 1.864266, -0.6305845,
    8.1104480e-2)                                 # polynomial for ξ
const impcfn_b = (0.2551543, -0.5455462, 0.3096816, 4.2568920e-2,
    -2.0123060e-2, -4.9607030e-3)                 # polynomial for φ
const impcfn_asymptotic_x = 2.0                   # asymptotic forms above this
const impcfn_asymptotic = (0.25, 1/32, 3/32)      # 1 + c₁/x ± c₂/x², c₃ = 3/32
const impcfn_small_x = 0.05                       # limiting forms below this
const impcfn_small_slope = 0.01917
const impcfn_two_exp_gamma = 1.1229               # 2 e^{-γ}

# Seaton's impact-parameter excitation rate (impactn)
const impactn_xm_max = 60.0
const impactn_psi_coeff = 1.644e5
const impactn_cr_coeff = 6.900e-5
const impactn_b0 = 10.0                           # start of the integration in b
const impactn_b_fraction = 100.0                  # step is b/this
const impactn_steps = 90
const impactn_wi_max = 100.0                      # stop when w/kT is this large
const impactn_cr_floor = 1e-20
const impactn_rel_tol = 1e-5
const impactn_inc_floor = 1e-7

# Simpson & Zhang (1988) semi-empirical excitation rates (szcoll)
const sz_abethe = (1.30, 0.59, 0.38, 0.286, 0.229, 0.192, 0.164, 0.141, 0.121,
    0.105, 0.100)
const sz_hbethe = (1.48, 3.64, 5.93, 8.32, 10.75, 12.90, 15.05, 17.20, 19.35,
    21.50, 2.15)
const sz_rbethe = (1.83, 1.60, 1.53, 1.495, 1.475, 1.46, 1.45, 1.45, 1.46,
    1.47, 1.48)
const sz_fvg1 = (1.133, 1.0785, 0.9935, 0.2328, -0.1296)
const sz_fvg2 = (-0.4059, -0.2319, 0.6282, -0.5598, 0.5299)
const sz_fvg3 = (0.07014, 0.02947, 0.3887, -1.181, 1.47)
const sz_fnn_coeff = 1.9603
const sz_cnn_coeff = 1.12
const sz_cnn_exp = 0.006

# Johnson (1972) hydrogenic excitation (erc)
const erc_y_max = 40.0
const erc_szcoll_ic = 10                          # szcoll for ionic charges from here
const erc_f_coeff = -1.2456e-10
const erc_z_coeff = 1.94
const erc_z_exp = 0.43
const erc_z_n1 = 0.45                             # n = 1
const erc_bn = (4.0, -18.63, 36.24, -28.09)
const erc_bn_n1 = -0.603                          # n = 1
const erc_s_coeff = 1.095e-10

# ln(n!) as a sum of logs (0 for n ≤ 0), as dfact()
logfact(n) = sum(log(i) for i in 1:n; init=0.0)

# hypergeometric function of integer parameters (hgf)
function hgf(ia, ib, ic, x)
    ser = 1.0
    hyp = 1.0
    for n in 0:min(-ia, -ib)
        ser = ser*(ia + n)*(ib + n)*x/((n + 1.0)*(ic + n))
        hyp += ser
    end
    hyp
end

"""
    expint(x)

`x e^x E₁(x)` (XSTAR's expint output `em1`): the Abramowitz & Stegun
polynomial fits for x > 1 and the series for x ≤ 1.
"""
function expint(x)
    if x <= 1
        a0, a1, a2, a3, a4, a5 = expint_series
        e1 = x > 0 ?
            a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5 - log(x) :
            -a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5 - log(-x)
        e1*x*expo(x)
    else
        b1, b2, b3, b4 = expint_denominator
        c1, c2, c3, c4 = expint_numerator
        (x^4 + c1*x^3 + c2*x^2 + c3*x + c4)/(x^4 + b1*x^3 + b2*x^2 + b3*x + b4)
    end
end

# first exponential integral E₁(t), as eint()
eint1(t) = expint(t)/t/expo(t)

"""
    anl1(ni, nf, lf, iq)

Spontaneous dipole rates `(alm, alp)` (s⁻¹) from `(ni, lf-1)` and `(ni, lf+1)` to
`(nf, lf)` for an ion of charge `iq`, from Gordon's formula.
"""
function anl1(ni, nf, lf, iq)
    K = constants()
    alm = 0.0
    alp = 0.0
    for li in (lf - 1, lf + 1)
        li < 0 && continue
        n, np, l = lf > li ? (nf, ni, lf) : (ni, nf, li)
        x1 = logfact(n + l)
        x2 = logfact(np + l - 1)
        x3 = logfact(2l - 1)
        x4 = logfact(n - l - 1)
        x5 = logfact(np - l)
        ia1 = -n + l + 1
        ia2 = ia1 - 2
        ib = -np + l
        ic = 2l
        x = -4.0*n*np/((n - np)*(n - np))
        y1 = hgf(ia1, ib, ic, x)
        y2 = hgf(ia2, ib, ic, x)
        rev = abs(n - np)
        rn = float(n + np)
        t = (l + 1)*log(4.0*n*np) + (rn - 2l - 2)*log(rev)
        t = t - log(4.0) - rn*log(rn)
        y1 = log(abs(y1 - y2*(rev/rn)^2)) + t
        t = expo(2y1 + x1 + x2 - 2x3 - x4 - x5)
        an = K.gordon_A*iq^4*max(li, lf)*t/(2.0*li + 1)
        an *= (1.0/nf^2 - 1.0/ni^2)^3
        li < lf && (alm = an)
        li > lf && (alp = an)
    end
    (alm, alp)
end

"""
    velimp(n, l, T, ic, z1, rm, ne, sum)

Rate (s⁻¹) for the l-changing collision `nl → nl-1` of a hydrogenic atom of
ionic charge `ic` with ions of charge `z1` and mass `rm` (electron masses), after
Pengelly & Seaton (1964), for temperature `T` (K), electron density `ne` and
total spontaneous decay rate `sum` of the level.
"""
function velimp(n, l, T, ic, z1, rm, ne, sum)
    K = constants()
    (l == 0 || sum == 0) && return 0.0
    den = l*(n^2 - l^2) + (l + 1)*(n^2 - (l + 1)^2)
    dnl = ps_dnl_coeff*z1/ic*z1/ic*n*n*(n*n - l*l - l - 1)
    pa = ps_rho_coeff/sum
    pd = K.ps_pd*sqrt(T/ne)
    alfa = K.ps_alfa*rm/T
    b = ps_b_coeff*sqrt(dnl)
    bb = b*b
    va = pd/pa
    vd = b/pd
    vb = sqrt(va*vd)
    ava, avb, avd = alfa*va^2, alfa*vb^2, alfa*vd^2
    xa, xb, xd = expo(-ava), expo(-avb), expo(-avd)
    ea = ava < ps_exp_cutoff ? expint(ava)/ava*xa : 0.0
    eb = expint(avb)/avb*xb
    ed = avd < ps_exp_cutoff ? expint(avd)/avd*xd : 0.0
    s = sqrt(π*alfa)
    if va > vd
        cn = avb > ps_series_cutoff ?
            s*(pa*pa*(2/alfa/alfa - xb*(vb^4 + 2vb*vb/alfa + 2/alfa/alfa)) +
                bb*xb + 2bb*eb - bb*ea) :
            s*bb*(1 + avb*(ps_series[1] - avb*ps_series[2]) + 2eb - ea)
    else
        ca = ava > ps_series_cutoff ?
            s*pa*pa*(2/alfa/alfa - xa*(va^4 + 2va*va/alfa + 2/alfa/alfa)) :
            s*pd*pd*va^4*alfa*(ps_series[1] - ava*ps_series[2] + ava*ava*ps_series[3])
        cad = s*pd*pd/alfa*(xa*(1 + ava) - xd*(1 + avd))
        cd = s*bb*(xd + ed)
        cn = ca + cad + cd
    end
    cn*l*(n^2 - l^2)/den
end

# Seaton's impact-parameter functions (impcfn): returns (xsi, phi)
function impcfn(x)
    a, b = impcfn_a, impcfn_b
    c1, c2, c3 = impcfn_asymptotic
    if x > impcfn_asymptotic_x
        xsi = π*x*exp(-2x)*(1 + c1/x + c2/x/x)
        phi = π/2*exp(-2x)*(1 + c1/x - c3/x/x)
    else
        xsi = 0.0
        phi = 0.0
        for n in 1:6
            xsi += a[n]*x^(n - 1)
            y = log(x)
            phi += b[n]*y^(n - 1)
        end
        x == 1 && (phi = b[1])
        if x < impcfn_small_x
            xsi = 1.0 + impcfn_small_slope/impcfn_small_x*x
            y = log(impcfn_two_exp_gamma/x)
            phi = y + x*x/4*(1 - 2y*y)
        end
    end
    (xsi, phi)
end

"""
    impactn(n, m, T, ic, amn)

Symmetric collision quantity `cmm` for the electron-impact transition `n → m` of
a hydrogenic ion of charge `ic`, by Seaton's impact-parameter method with the
strong-coupling cross sections only.
"""
function impactn(n, m, T, ic, amn)
    K = constants()
    xm = K.Ry_K*ic*ic/T/m/m
    xm > impactn_xm_max && return 0.0
    rm, z1 = 1.0, 1.0
    tk = K.kB_eV*T
    ecm = K.Rinf_invcm*ic*ic*(1.0/n/n - 1.0/m/m)
    ecm3 = ecm^3
    ecm = -ecm
    psi = impactn_psi_coeff*amn/ecm3
    cr, fi, wo, b = 0.0, 0.0, 0.0, impactn_b0
    ev = abs(ecm)/K.invcm_per_eV
    done = false
    while !done
        del = b/impactn_b_fraction
        for _ in 1:impactn_steps
            b -= del
            xsi, phi = impcfn(b)
            w = ic*rm*ev/b*sqrt(2*xsi*psi)
            wi = w + ecm/K.invcm_per_eV/2
            if wi/tk >= impactn_wi_max
                done = true
                break
            end
            wi <= 0 && continue
            ff = (xsi/2 + phi)*exp(-wi/tk)
            crinc = (fi + ff)/2*(wi - wo)
            cr += crinc
            cr < impactn_cr_floor && continue
            fi, wo = ff, wi
            if crinc/cr < impactn_rel_tol && crinc > impactn_inc_floor
                done = true
                break
            end
        end
    end
    cr = impactn_cr_coeff*z1*z1*sqrt(rm/T)*psi*cr/tk
    cr*m*m*exp(xm)
end

# Simpson & Zhang (1988) semi-empirical excitation rate, n=ni → nj (szcoll)
function szcoll(ni, nj, T, ic)
    K = constants()
    rn2 = (float(ni)/float(nj))^2
    g1 = g2 = g3 = 0.0
    if ni == 1
        g1, g2, g3 = sz_fvg1[1], sz_fvg2[1], sz_fvg3[1]
    elseif ni == 2
        g1, g2, g3 = sz_fvg1[2], sz_fvg2[2], sz_fvg3[2]
    else
        g1 = sz_fvg1[3] + sz_fvg1[4]/ni + sz_fvg1[5]/ni/ni
        g2 = -(sz_fvg2[3] + sz_fvg2[4]/ni + sz_fvg2[5]/ni/ni)/ni
        g3 = (sz_fvg3[3] + sz_fvg3[4]/ni + sz_fvg3[5]/ni/ni)/ni/ni
    end
    xx = 1 - rn2
    gaunt = g1 + g2/xx + g3/xx/xx
    fnn = sz_fnn_coeff*gaunt/xx^3*ni/nj^3
    if ni < length(sz_abethe)
        an, hn, rrn = sz_abethe[ni], sz_hbethe[ni], sz_rbethe[ni]
    else
        an, hn, rrn = sz_abethe[end]/ni, sz_hbethe[end]*ni, sz_rbethe[end]
    end
    ann = fnn*4*ni^4/(1 - rn2)
    dnn = ann*hn*(xx^rrn - an*rn2)
    cnn = sz_cnn_coeff*ni*ann*xx
    nj - ni == 1 && (cnn *= expo(-sz_cnn_exp*(ni - 1)^6/ic))
    yy = K.Ry_K_sz*ic*ic*(1/float(ni*ni) - 1/float(nj*nj))/T
    e1 = eint1(yy)
    K.sz_rate_coeff/sqrt(T)/ni/ni/ic/ic*(dnn*expo(-yy) + (ann + yy*(cnn - dnn))*e1)
end

"""
    erc(n, m, T, ic, a)

Hydrogenic electron-impact excitation and de-excitation rates `(se, sd)`
(cm³ s⁻¹) between principal quantum numbers `n < m` for an ion of charge `ic`
at `T` (K); `a` is the summed spontaneous rate used by the impact-parameter fit.
"""
function erc(n, m, T, ic, a)
    erc_ryd_K = constants().Ry_K_erc
    if ic != 1
        ym = erc_ryd_K*ic*ic/T/m/m
        if ic < erc_szcoll_ic
            ym > erc_y_max && return (0.0, 0.0)
            sm = impactn(n, m, T, ic, a)
            xn = 1.0/n/n - 1.0/m/m
            yn = erc_ryd_K*ic*ic*xn/T
            s = sm/n/n/expo(ym)
            sd = s*n/m*n/m
            se = yn < erc_y_max ? s*expo(-yn) : 0.0
        else
            se = szcoll(n, m, T, ic)
            xn = 1.0/n/n - 1.0/m/m
            yn = erc_ryd_K*ic*ic*xn/T
            sd = se*n*n/m/m*expo(yn)
        end
        (se, sd)
    else
        xn = 1.0/n/n - 1.0/m/m
        f = erc_f_coeff*a/xn/xn
        yn = erc_ryd_K*xn/T
        ym = erc_ryd_K/T/m/m
        z = n == 1 ? yn + erc_z_n1*xn : erc_z_coeff*xn*n^erc_z_exp + yn
        dif = z - yn
        e1y = expint(yn)
        e1z = expint(z)
        e2 = (1 - e1y)/yn - expo(-dif)*(1 - e1z)/z
        rn, rm = float(n), float(m)
        ann = -2*f*m*m/xn/n/n
        bn = n == 1 ? erc_bn_n1 :
            (erc_bn[1] + erc_bn[2]/n + erc_bn[3]/n/n + erc_bn[4]/n^3)/n
        bnn = (1 + 4/(xn*n*n*3) + bn/(rn^4*xn*xn))*4/(rm^3*xn*xn)
        s = ann*((1/yn + 0.5)*e1y/yn - (1/z + 0.5)*e1z*expo(-dif)/z)
        s += e2*(bnn - ann*log(2/xn))
        s = erc_s_coeff*yn*yn*sqrt(T)*s/xn
        (s*expo(-yn), s*n/m*n/m)
    end
end

# ---------------------------------------------------------------------------
# collisional ionization of hydrogenic ions: eint, szirc and irc

const eint_floor = 1e-34                          # eint divides by at least this

"""
    eint(t)

The exponential integrals `(E₁, E₂, E₃)(t)` as XSTAR's `eint`.
"""
function eint(t)
    e1 = expint(t)/max(eint_floor, t*expo(t))
    e2 = exp(-t) - t*e1
    e3 = (expo(-t) - t*e2)/2
    (e1, e2, e3)
end

# szirc: Bethe-type fit parameters for n = 1:10 and the scalings for higher n (single-precision literals)
const szirc_a = (1.134, 0.603, 0.412, 0.313, 0.252, 0.211, 0.181, 0.159, 0.142, 0.128, 1.307)
const szirc_h = (1.48, 3.64, 5.93, 8.32, 10.75, 12.90, 15.05, 17.20, 19.35, 21.50, 2.15)
const szirc_r = (2.20, 1.90, 1.73, 1.65, 1.60, 1.56, 1.54, 1.52, 1.52, 1.52, 1.52)
const szirc_nmax = 11
const szirc_coeff = Float64(4.6513f-3)
const szirc_exchange = 3.36

"""
    szirc(n, T, rz, rno)

Collisional ionization coefficient from level `n` for the effective charge `rz` and highest bound level
`rno` at the temperature `T` (K), by Sampson and Zhang's fit (XSTAR's `szirc`).
"""
function szirc(n, T, rz, rno)
    rc = Float64(trunc(Int, rno))
    an, hn, rrn = n < szirc_nmax ?
        (szirc_a[n], szirc_h[n], szirc_r[n]) :
        (szirc_a[end]/n, szirc_h[end]*n, szirc_r[end])
    K = constants()
    tt = T*K.kB_szirc
    rn = Float64(n)
    yy = rz*rz*K.Ry_erg_szirc/tt*(1/rn/rn - 1/rc/rc - (1/(rc - 1)^2 - 1/rc/rc)/4)
    e1, e2, e3 = eint(yy)
    szirc_coeff*sqrt(tt)*rn^5/rz^4*an*yy*(
        e1/rn - (exp(-yy) - yy*e3)/(3rn) +
        (yy*e2 - 2yy*e1 + exp(-yy))*3hn/rn/(3 - rrn) + (e1 - e2)*szirc_exchange*yy)
end

# irc: Hydrogen-like fits for the ionization from n = 1, 2, and higher (rc = 1)
const irc_fit1 = ((1.133, 0.4059, 0.07014), (1.0785, 0.2319, 0.02947))
const irc_rn = (0.45, 0.653)
const irc_coeff = 1.095e-10

"""
    irc(n, T, rc, rno)

Collisional ionization coefficient from level `n` of an ion of effective charge `rc` for the
highest bound level `rno` at the temperature `T` (K) (XSTAR's `irc`). For `rc ≠ 1` it is `szirc`.
"""
function irc(n, T, rc, rno)
    rc != 1 && return szirc(n, T, rc, rno)
    xo = 1 - n*n/rno/rno
    yn = xo*constants().Ry_K_erc/(T*n*n)
    if n < 3
        c = irc_fit1[n]
        an = 1.9603*n*(c[1]/3/xo^3 - c[2]/4/xo^4 + c[3]/5/xo^5)
        bn = n == 1 ? 2/3*n*n/xo*(3 + 2/xo - 0.603/xo/xo) :
            2/3*n*n/xo*(3 + 2/xo + (4 - 18.63/n + 36.24/(n*n) - 28.09/(n*n*n))/n/xo/xo)
        rn = irc_rn[n]
    else
        g0 = (0.9935 + 0.2328/n - 0.1296/(n*n))/3/xo^3
        g1 = -(0.6282 - 0.5598/n + 0.5299/(n*n))/(n*4)/xo^4
        g2 = (0.3887 - 1.181/n + 1.470/(n*n))/(n*n*5)/xo^5
        an = 1.9603*n*(g0 + g1 + g2)
        bn = (4 - 18.63/n + 36.24/(n*n) - 28.09/(n*n*n))/n
        bn = (3 + 2/xo + bn/xo/xo)*2*n*n/3/xo
        rn = 1.94*n^(-1.57)
    end
    rn *= xo
    zn = rn + yn
    ey, ez = expint(yn), expint(zn)
    se = an*(ey/yn/yn - exp(-rn)*ez/zn/zn)
    ey = 1 + 1/yn - ey*(2/yn + 1)
    ez = exp(-rn)*(1 + 1/zn - ez*(2/zn + 1))
    se += (bn - an*log(2*n*n/xo))*(ey - ez)
    se*sqrt(T)*yn*yn*n*n*irc_coeff/xo
end

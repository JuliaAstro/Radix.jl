# Hydrogenic collision routines used by XSTAR's ucalc type 63 (CollisionProb):
# Gordon's radial integrals (anl1), Pengelly & Seaton / Hummer & Storey l-changing
# rates (velimp), and Johnson / Seaton n-changing excitation (erc, impactn,
# szcoll). Translated from xstarsub.f; the arithmetic is Float64 where the Fortran
# is single precision.

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
        a0, a1, a2, a3, a4, a5 = -0.57721566, 0.99999193, -0.24991055,
            0.05519968, -0.00976004, 0.00107857
        e1 = x > 0 ?
            a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5 - log(x) :
            -a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5 - log(-x)
        e1*x*expo(x)
    else
        b1, b2, b3, b4 = 9.5733223454, 25.6329561486, 21.0996530827, 3.9584969228
        c1, c2, c3, c4 = 8.5733287401, 18.0590169730, 8.6347608925, 0.2677737343
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
        an = 2.6761e9*iq^4*max(li, lf)*t/(2.0*li + 1)
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
    (l == 0 || sum == 0) && return 0.0
    den = l*(n^2 - l^2) + (l + 1)*(n^2 - (l + 1)^2)
    dnl = 6.0*z1/ic*z1/ic*n*n*(n*n - l*l - l - 1)
    pa = 0.72/sum
    pd = 6.90*sqrt(T/ne)
    alfa = 3.297e-12*rm/T
    b = 1.157*sqrt(dnl)
    bb = b*b
    va = pd/pa
    vd = b/pd
    vb = sqrt(va*vd)
    ava, avb, avd = alfa*va^2, alfa*vb^2, alfa*vd^2
    xa, xb, xd = expo(-ava), expo(-avb), expo(-avd)
    ea = ava < 50 ? expint(ava)/ava*xa : 0.0
    eb = expint(avb)/avb*xb
    ed = avd < 50 ? expint(avd)/avd*xd : 0.0
    s = sqrt(π*alfa)
    if va > vd
        cn = avb > 1e-3 ?
            s*(pa*pa*(2/alfa/alfa - xb*(vb^4 + 2vb*vb/alfa + 2/alfa/alfa)) +
                bb*xb + 2bb*eb - bb*ea) :
            s*bb*(1 + avb*(1/3 - avb/4) + 2eb - ea)
    else
        ca = ava > 1e-3 ?
            s*pa*pa*(2/alfa/alfa - xa*(va^4 + 2va*va/alfa + 2/alfa/alfa)) :
            s*pd*pd*va^4*alfa*(1/3 - ava/4 + ava*ava/10)
        cad = s*pd*pd/alfa*(xa*(1 + ava) - xd*(1 + avd))
        cd = s*bb*(xd + ed)
        cn = ca + cad + cd
    end
    cn*l*(n^2 - l^2)/den
end

# Seaton's impact-parameter functions (impcfn): returns (xsi, phi)
function impcfn(x)
    a = (0.9947187, 0.6030883, -2.372843, 1.864266, -0.6305845, 8.1104480e-2)
    b = (0.2551543, -0.5455462, 0.3096816, 4.2568920e-2, -2.0123060e-2, -4.9607030e-3)
    if x > 2
        xsi = π*x*exp(-2x)*(1 + 0.25/x + 1/32/x/x)
        phi = π/2*exp(-2x)*(1 + 0.25/x - 3/32/x/x)
    else
        xsi = 0.0
        phi = 0.0
        for n in 1:6
            xsi += a[n]*x^(n - 1)
            y = log(x)
            phi += b[n]*y^(n - 1)
        end
        x == 1 && (phi = b[1])
        if x < 0.05
            xsi = 1.0 + 0.01917/0.05*x
            y = log(1.1229/x)
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
    xm = 157888.0*ic*ic/T/m/m
    xm > 60 && return 0.0
    rm, z1 = 1.0, 1.0
    tk = 8.617e-5*T
    ecm = 109737.0*ic*ic*(1.0/n/n - 1.0/m/m)
    ecm3 = ecm^3
    ecm = -ecm
    psi = 1.644e5*amn/ecm3
    cr, fi, wo, b = 0.0, 0.0, 0.0, 10.0
    ev = abs(ecm)/8065.48
    done = false
    while !done
        del = b/100
        for _ in 1:90
            b -= del
            xsi, phi = impcfn(b)
            w = ic*rm*ev/b*sqrt(2*xsi*psi)
            wi = w + ecm/8065.48/2
            if wi/tk >= 100
                done = true
                break
            end
            wi <= 0 && continue
            ff = (xsi/2 + phi)*exp(-wi/tk)
            crinc = (fi + ff)/2*(wi - wo)
            cr += crinc
            cr < 1e-20 && continue
            fi, wo = ff, wi
            if crinc/cr < 1e-5 && crinc > 1e-7
                done = true
                break
            end
        end
    end
    cr = 6.900e-5*z1*z1*sqrt(rm/T)*psi*cr/tk
    cr*m*m*exp(xm)
end

# Simpson & Zhang (1988) semi-empirical excitation rate, n=ni → nj (szcoll)
function szcoll(ni, nj, T, ic)
    abethe = (1.30, 0.59, 0.38, 0.286, 0.229, 0.192, 0.164, 0.141, 0.121, 0.105, 0.100)
    hbethe = (1.48, 3.64, 5.93, 8.32, 10.75, 12.90, 15.05, 17.20, 19.35, 21.50, 2.15)
    rbethe = (1.83, 1.60, 1.53, 1.495, 1.475, 1.46, 1.45, 1.45, 1.46, 1.47, 1.48)
    fvg1 = (1.133, 1.0785, 0.9935, 0.2328, -0.1296)
    fvg2 = (-0.4059, -0.2319, 0.6282, -0.5598, 0.5299)
    fvg3 = (0.07014, 0.02947, 0.3887, -1.181, 1.47)
    eion, c0 = 1.578203e5, 8.63e-6
    rn2 = (float(ni)/float(nj))^2
    delnn = (1/float(ni*ni) - 1/float(nj*nj))^2
    g1 = g2 = g3 = 0.0
    if ni == 1
        g1, g2, g3 = fvg1[1], fvg2[1], fvg3[1]
    elseif ni == 2
        g1, g2, g3 = fvg1[2], fvg2[2], fvg3[2]
    else
        g1 = fvg1[3] + fvg1[4]/ni + fvg1[5]/ni/ni
        g2 = (fvg2[3] + fvg2[4]/ni + fvg2[5]/ni/ni)/ni*(-1.0)
        g3 = (fvg3[3] + fvg3[4]/ni + fvg3[5]/ni/ni)/ni/ni
    end
    xx = 1 - rn2
    gaunt = g1 + g2/xx + g3/xx/xx
    fnn = 1.9603*gaunt/xx^3*ni/nj^3
    if ni < 11
        an, hn, rrn = abethe[ni], hbethe[ni], rbethe[ni]
    else
        an, hn, rrn = abethe[11]/ni, hbethe[11]*ni, rbethe[11]
    end
    ann = fnn*4*ni^4/(1 - rn2)
    dnn = ann*hn*(xx^rrn - an*rn2)
    cnn = 1.12*ni*ann*xx
    nj - ni == 1 && (cnn *= expo(-0.006*(ni - 1)^6/ic))
    yy = eion*ic*ic*(1/float(ni*ni) - 1/float(nj*nj))/T
    e1 = eint1(yy)
    c0/sqrt(T)/ni/ni/ic/ic*(dnn*expo(-yy) + (ann + yy*(cnn - dnn))*e1)
end

"""
    erc(n, m, T, ic, a)

Hydrogenic electron-impact excitation and de-excitation rates `(se, sd)`
(cm³ s⁻¹) between principal quantum numbers `n < m` for an ion of charge `ic`
at `T` (K); `a` is the summed spontaneous rate used by the impact-parameter fit.
"""
function erc(n, m, T, ic, a)
    if ic != 1
        ym = 157803.0*ic*ic/T/m/m
        if ic < 10
            ym > 40 && return (0.0, 0.0)
            sm = impactn(n, m, T, ic, a)
            xn = 1.0/n/n - 1.0/m/m
            yn = 157803.0*ic*ic*xn/T
            s = sm/n/n/expo(ym)
            sd = s*n/m*n/m
            se = yn < 40 ? s*expo(-yn) : 0.0
        else
            se = szcoll(n, m, T, ic)
            xn = 1.0/n/n - 1.0/m/m
            yn = 157803.0*ic*ic*xn/T
            sd = se*n*n/m/m*expo(yn)
        end
        (se, sd)
    else
        xn = 1.0/n/n - 1.0/m/m
        f = -1.2456e-10*a/xn/xn
        yn = 157803.0*xn/T
        ym = 157803.0/T/m/m
        z = n == 1 ? yn + 0.45*xn : 1.94*xn*n^0.43 + yn
        dif = z - yn
        e1y = expint(yn)
        e1z = expint(z)
        e2 = (1 - e1y)/yn - expo(-dif)*(1 - e1z)/z
        rn, rm = float(n), float(m)
        ann = -2*f*m*m/xn/n/n
        bn = n == 1 ? -0.603 : (4 - 18.63/n + 36.24/n/n - 28.09/n^3)/n
        bnn = (1 + 4/(xn*n*n*3) + bn/(rn^4*xn*xn))*4/(rm^3*xn*xn)
        s = ann*((1/yn + 0.5)*e1y/yn - (1/z + 0.5)*e1z*expo(-dif)/z)
        s += e2*(bnn - ann*log(2/xn))
        s = 1.095e-10*yn*yn*sqrt(T)*s/xn
        (s*expo(-yn), s*n/m*n/m)
    end
end

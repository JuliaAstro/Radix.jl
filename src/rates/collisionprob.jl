# Collisional transition probability Cik for N-electron ion computed by
# quantum defect theory (or hydrogenic): i1= 1; i2= i (lower level);
# i3= k (upper level); i4= Z; i5= ionN

# XSTAR data type: 63

const CollisionProbDesc = "h-like cij, bautista (hlike ion)"

struct CollisionProb{I} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    transition::Transition{I}   # the two levels, in stored order
    Z::I      # atomic number
    ion::I    # ion index (XSTAR ionN)
end

function CollisionProb(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    CollisionProb(Int8(rate), label, Transition(ivec[2], ivec[3]), ivec[4:5]...)
end

"""
    rate(coef::CollisionProb, cell; levels, index=false)

Collisional transition rates between two levels of a hydrogenic ion (XSTAR
ucalc type 63), for `nl → n'l'` with `|Δl| = 1`. `levels` is a `level_table`.
`init` and `final` are the two levels in stored order, `frate` is the rate
`init → final` and `irate` the reverse (s⁻¹, including the factor `nₑ`). Both
are 0 when the levels do not have `|Δl| = 1`, when `ΔE/kT > 50`, or if a level
is missing from `levels`.

For `n = n'` the rate is the l-changing rate of Pengelly & Seaton (`velimp`) with
the detailed-balance partner. For `n ≠ n'` it is the electron-impact rate `erc`
distributed over `l` with Gordon's radiative rates (`anl1`). XSTAR applies the
large (de-excitation) rate of this branch to the *upward* direction (the source
comments "check if ans1 and ans2 are correct or inverted"); this is reproduced
as is.
"""
function rate(coef::CollisionProb, cell::Cell; levels, index=false,
    verbose=false)

    none = (; init=0, final=0, frate=0., irate=0.)
    a = get(levels, (coef.ion, coef.transition.lower), nothing)
    b = get(levels, (coef.ion, coef.transition.upper), nothing)
    (a === nothing || b === nothing) && return none
    elin = 12398.54/abs(Float64(b.E) - Float64(a.E) + 1e-24)
    12398.54/elin/(0.861707*cell.T) > 50 && return none
    index && return (; none..., init=a.level, final=b.level)

    ni, li, nf, lf, Z = Int(a.n), Int(a.L), Int(b.n), Int(b.L), Int(coef.Z)
    T = cell.T*1e4                                    # K
    ans1 = ans2 = 0.0
    if abs(lf - li) == 1
        if nf == ni
            lii = max(lf, li)
            sum = 0.0
            for nn in max(1, min(lf, li)):ni - 1
                lii >= 1 && (sum += anl1(ni, nn, lii - 1, Z)[2])
                nn > lii + 1 && (sum += anl1(ni, nn, lii + 1, Z)[1])
            end
            # XSTAR's amcrs is called with ecm = 0, which always ends in velimp
            cn = velimp(ni, lii, T, Z, 1.0, 1800.0, cell.nₑ, sum)
            if lf < li
                ans1, ans2 = cn, cn*a.g/b.g
            else
                ans1, ans2 = cn*b.g/a.g, cn
            end
        else
            nu, nll = max(ni, nf), min(ni, nf)
            sum = 0.0
            aa1 = 0.0
            for lff in 0:nll - 1
                alm, alp = anl1(nu, nll, lff, Z)
                sum += alp*(2lff + 3)
                lff > 0 && (sum += alm*(2lff - 1))
                lff == lf && li > lf && (aa1 = alp)
                lff == lf && li < lf && (aa1 = alm)
            end
            se, sd = erc(nll, nu, T, Z, sum)
            ans1 = se*(2lf + 1)*aa1/sum
            ans2 = sd*(2li + 1)*aa1/sum
            (nf > ni || lf > li) && ((ans1, ans2) = (ans2, ans1))
        end
    end
    (; init=a.level, final=b.level, frate=ans1*cell.nₑ, irate=ans2*cell.nₑ)
end

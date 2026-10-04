# Coefficients for recombination and photoionization cross sections of
# superlevels: r1−rjnd = (ne( j), j= 1, jnd); rjnd+1−rjnd+nt = (Te( j),
# j= 1, jnt); rjnd+nt+1−rjnd+nt+nt*nd = ((log α( j, j′), j′ = 1, j′ nd),
# j= 1, jnt); rjnd+nt+nt*nd+1−rjnd+nt+nt*nd+2*nx = (E( j),σ( j), j= 1, jnx);
# i1= nd; i2= nt; i3= nx; i4= n; i5= L; i6= 2S + 1; i7= Z; i8= kN−1;
# i9= ionN−1; i10= iN ; i11= ionN

# XSTAR data type: 70

const PhotoionizeSuperDesc = "Coefficients for phot x-section of suplevels"

struct PhotoionizeSuper{I, R, L<:LevelTable} <: AbstractRate
    rtype::Int8                 # XSTAR rate type (lrtyp)
    label::String
    n::I                # principal quantum number
    L::I                # orbital angular momentum
    spin_mult::I        # 2S+1
    Z::I                # atomic number
    parent::Parent{I}           # parent ion and level
    level::I            # level index
    ion::I              # ion index (XSTAR ionN)
    ne_grid::Vector{R}  # log₁₀ densities (cm⁻³)
    T_grid::Vector{R}   # log₁₀ temperatures (K)
    logα::Matrix{R}     # log₁₀ α (cm³ s⁻¹), size (T, ne)
    E_grid::Vector{R}   # energies
    σ::Vector{R}        # cross sections
    levels::L              # the level data of the database (a LevelTable)
end

function PhotoionizeSuper(rate::Int32, label::String, ivec::I, rvec::R, levels::LevelTable) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    Nd, Nt, Nx = ivec[1:3]
    # ucalc reads the table with the temperature running fastest
    PhotoionizeSuper(Int8(rate), label, ivec[4:7]..., Parent(ivec[9], ivec[8]), ivec[10:11]..., Vector(rvec[1:Nd]),
        Vector(rvec[Nd+1:Nd+Nt]), reshape(rvec[Nd+Nt+1:Nd+Nt+Nd*Nt], Int(Nt), Int(Nd)),
        Vector(rvec[Nd+Nt+Nd*Nt+1:2:end-1]), Vector(rvec[Nd+Nt+Nd*Nt+2:2:end]), levels)
end

"""
    superlevel_cross_section(coef::PhotoionizeSuper, T, den, E_th)

Recombination coefficient `rec` and the cross-section table scaled to it (XSTAR's `calt70`): `rec` is
interpolated in log₁₀ T and log₁₀ n (linearly in T, quadratically in n between three densities), the
table is scaled so that the Milne relation at `T` (K) with the threshold `E_th` (Ry) gives `rec`, capped
at 10⁶ Mb, and its tail below 10⁻⁶ of the first point is dropped. As in ucalc, the second tabulated
log₁₀ density is limited to 8 and the density to the last tabulated one. Returns `rec`, `E_grid` and `σ`.
"""
function superlevel_cross_section(coef::PhotoionizeSuper, T, den, E_th)
    logne = Float64.(coef.ne_grid)
    logT = Float64.(coef.T_grid)
    logα = Float64.(coef.logα)
    nden, ntem = length(logne), length(logT)
    nden > 1 && (logne[2] = min(logne[2], super_log_density_max))
    rne, rte = log10(den), log10(T)

    in = 1                                          # density interval
    if nden > 1
        rne = min(rne, logne[nden])
        if rne > logne[1]
            in = trunc(Int, rne/logne[nden]*nden) - 1
            in >= nden && (in -= 1)
            while true
                in += 1
                in < nden && rne >= logne[in + 1] && continue
                rne < logne[in] && (in -= 2; continue)
                break
            end
        end
    end
    it = 1                                          # temperature interval
    if rte >= logT[1]
        dt = (logT[ntem] - logT[1])/ntem
        it = trunc(Int, (rte - logT[1])/dt)
        while true
            it += 1
            if it >= ntem
                it = ntem - 1
            elseif rte >= logT[it + 1]
                continue
            elseif rte < logT[it]
                it -= 2
                continue
            end
            break
        end
    end

    at(j) = logα[it, j] + (logα[it + 1, j] - logα[it, j])/(logT[it + 1] - logT[it])*(rte - logT[it])
    rec = at(in)
    if nden > 2 && 1 < in < nden
        x1, x2, x3 = logne[in - 1], logne[in], logne[in + 1]
        y1, y2, y3 = at(in - 1), rec, at(in + 1)
        denom = (x1*x1 - x2*x2)*(x1 - x3) - (x1*x1 - x3*x3)*(x1 - x2)
        a = ((y1 - y2)*(x1 - x3) - (y1 - y3)*(x1 - x2))/denom
        b = -((y1 - y2)*(x1*x1 - x3*x3) - (y1 - y3)*(x1*x1 - x2*x2))/denom
        c = y1 - a*x1*x1 - b*x1
        rec = a*rne*rne + b*rne + c
    elseif nden > 1 && in < nden                    # (at the last density rne equals logne[in])
        rec += (at(in + 1) - rec)/(logne[in + 1] - logne[in])*(rne - logne[in])
    end
    rec = exp10(rec)

    ε = Float64.(coef.E_grid)
    xs = Float64.(coef.σ)
    scale = rec/(milne_flat + milne_recombination(T, ε, xs, E_th))
    xs = min.(xs.*scale, super_xs_cap)
    last = length(xs)
    for i in eachindex(xs)
        xs[i] > xs[1]*super_xs_trim && (last = i)
    end
    (rec, ε[1:last], xs[1:last])
end

"""
    rate(coef::PhotoionizeSuper, cell; neutral=false, index=false)

Photoionization and recombination of a superlevel from the tabulated recombination coefficients and
cross section (XSTAR ucalc type 70, "old type 70"), see `photoionize_superlevel`. The level data comes from the coefficient (its `levels` field) and `radiation` is the spectrum (none by default). `neutral=true` limits the density to
10⁸ cm⁻³ in the interpolation, as ucalc does for the first ion of an element. The other keywords of the
photoionization rates are accepted and ignored.
"""
function rate(coef::PhotoionizeSuper, cell::Cell; radiation=NO_RADIATION, neutral=false, index=false, ptmp=nothing, abund=nothing, lfast=nothing, opacity=nothing)
    photoionize_superlevel(coef, cell,
        (T, n, E_th) -> superlevel_cross_section(coef, T, neutral ? min(n, super_density_cap) : n, E_th);
        radiation, index)
end

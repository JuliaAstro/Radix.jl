# Coefficients for recombination and photoionization cross sections of
# superlevels: r1−rjnd = (ne( j), j=1, jnd); rjnd+1−rjnd+nt = (Te( j), j= 1,
# jnt); rjnd+nt+1−rjnd+nt+nt*nd = ((α( j, j′), j′= 1, j′nd), j= 1, jnt);
# rjnd+nt+nt*nd+1−rjnd+nt+nt*nd+2*nx = (E( j),σ( j), j= 1, jnx); i1= nd;
# i2= nt; i3= nx; i4= n; i5= L; i6= 2S + 1; i7= Z; i8= kN−1; i9= ionN−1;
# i10= iN ; i11= ionN

# XSTAR data type: 99

const PhotoRecombXDesc = "Recombination and photoionization of superlevels"

struct PhotoRecombX{I, R} <: AbstractRate
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
    α::Matrix{R}        # α (cm³ s⁻¹; values below -1e-31 are log₁₀ α), size (ne, T)
    E_grid::Vector{R}   # energies
    σ::Vector{R}        # cross sections
end

function PhotoRecombX(rate::Int32, label::String, ivec::I, rvec::R) where
    {I<:AbstractVector{Int32}, R<:AbstractVector{Float32}}

    nd, nt, nx = ivec[1:3]
    PhotoRecombX(Int8(rate), label, ivec[4:7]..., Parent(ivec[9], ivec[8]), ivec[10:11]..., Vector(rvec[1:nd]),
        Vector(rvec[nd+1:nd+nt]), reshape(rvec[nd+nt+1:nd+nt+nd*nt], Int(nd), Int(nt)),
        Vector(rvec[nd+nt+nd*nt+1:2:nd+nt+nd*nt+2*nx-1]),
        Vector(rvec[nd+nt+nd*nt+2:2:nd+nt+nd*nt+2*nx]))
end

const photorecombx_floor = 1f-30      # added to α before taking the logarithm
const photorecombx_log_min = -1f-31   # values below this are already logarithms
const photorecombx_T_inner = (0.999f0, 1.001f0)   # factors that keep log₁₀ T inside the table

"""
    superlevel_cross_section(coef::PhotoRecombX, T, den, E_th)

Recombination coefficient `rec` and the cross-section table scaled to it (XSTAR's `calt99`, in single
precision): log₁₀ α is interpolated linearly in log₁₀ T and in log₁₀ n between the two densities that
bracket `den`; at the lowest density, outside the table and at the highest density the lowest
density column is used. log₁₀ T is limited to the table. The table is scaled so that the Milne relation at
`T` (K) with the threshold `E_th` (Ry) gives `rec`. Returns `rec`, `E_grid` and `σ` (Mb).
"""
function superlevel_cross_section(coef::PhotoRecombX, T, den, E_th)
    dens, temps = Float32.(coef.ne_grid), Float32.(coef.T_grid)
    nden, ntem = length(dens), length(temps)
    logα = map(a -> a > photorecombx_log_min ? log10(a + photorecombx_floor) : a, Float32.(coef.α))   # (nden, ntem)
    rne, rte = log10(Float32(den)), log10(Float32(T))

    in = 1
    if nden > 1 && rne > dens[1]
        in = 0
        for i in 1:nden - 1
            rne >= dens[i] && rne <= dens[i + 1] && (in = i)
        end
    end
    in = max(in, 1)
    if rte < temps[1] || rte > temps[ntem]
        rte = min(photorecombx_T_inner[1]*temps[ntem], max(photorecombx_T_inner[2]*temps[1], rte))
    end
    it = 1
    for i in 1:ntem - 1
        rte >= temps[i] && rte < temps[i + 1] && (it = i)
    end

    at(j) = logα[j, it] + (logα[j, it + 1] - logα[j, it])/(temps[it + 1] - temps[it])*(rte - temps[it])
    rec1 = at(in)
    rr = rec1
    if in != nden && in > 1
        rr = rec1 + (at(in + 1) - rec1)*(1/(dens[in + 1] - dens[in]))*(rne - dens[in])
    end
    rec = exp10(rr)

    ε = Float64.(coef.E_grid)
    xs = Float64.(coef.σ)
    scale = rec/Float32(milne_recombination(T, ε, xs, E_th))
    (Float64(rec), ε, xs.*Float64(scale))
end

"""
    rate(coef::PhotoRecombX, cell; levels, radiation, nlev, index=false)

Photoionization and recombination of a superlevel from the tabulated recombination coefficients and
cross section (XSTAR ucalc type 99), see `photoionize_superlevel`. The energies measured from the levels
are corrected with the energy of the levels. The other keywords of the photoionization rates are
accepted and ignored.
"""
function rate(coef::PhotoRecombX, cell::Cell; levels, radiation, nlev, index=false,
    ptmp=nothing, abund=nothing, lfast=nothing, opacity=nothing, verbose=false)

    photoionize_superlevel(coef, cell; levels, radiation, nlev, index, correct_energy=true,
        tabulate=(T, n, E_th) -> superlevel_cross_section(coef, T, n, E_th))
end

# Sample real records of one XSTAR data type from atdb.fits and write the input
# for the ucalc driver (drvu): one case per record and (t, nH, xe) condition,
# with the level data of the record's ion.
#
#   julia generate.jl /path/to/atdb.fits OUTDIR TYPE [NRECORDS]
#
# then:  drvu < OUTDIR/typeNN.in > OUTDIR/typeNN.out
using FITSFiles, Random

const RADIATION_TYPES = (49, 50, 53, 59, 74, 85, 88) # types that need a spectrum
const PARENT_TYPES = (49, 53, 59)    # types that leave a level of the next ion
# positions of the parent level and the parent ion among the integers
parent_positions(type) = type == 59 ? (4, 5) : (5, 6)
const CFRAC = (0.0, 0.25, 0.0, 0.5)  # covering fraction per condition (line rates)
const MAX_TABLE = 800                # longest cross-section table (reals) sampled

const LFAST = (1, 2, 3, 3)             # ucalc's speed switch, per condition (radiation types)

const CONDITIONS = (            # t [1e4 K], xpx [cm⁻³], xee, xh0/xpx, xh1/xpx
    (0.3, 1e4, 0.5, 0.5, 0.5),
    (1.0, 1e8, 0.9, 0.1, 0.9),
    (10.0, 1e10, 1.0, 1e-3, 1.0),
    (100.0, 1e12, 1.0, 1e-6, 1.0),
)

# optional: TAG and SELECTION "ion:level:kpar,ion:level:kpar,..." to sample chosen records only
function main(path, outdir, type, nrec = 25, tag = "", selection = "")
    a = open(io -> fits(io; scale = false), path)
    ptr  = reshape(a[2].data[:pointers][1], 10, :)
    rdat = a[3].data[:reals][1]
    idat = a[4].data[:integers][1]
    N = findfirst(j -> ptr[2, j] == 0, 1:size(ptr, 2)) - 1
    ints(j)  = Int.(idat[ptr[9, j]:ptr[9, j] + ptr[6, j] - 1])
    reals(j) = Float64.(rdat[ptr[8, j]:ptr[8, j] + ptr[5, j] - 1])

    # levels per ion: ion => [(level, E, g, Einf, n, 2S+1, L)]
    levels = Dict{Int, Vector{NTuple{7, Float64}}}()
    for j in 1:N
        ptr[2, j] == 6 || continue
        i, r = ints(j), reals(j)                       # n, 2S+1, L, Z, level, ion ; E, g, ν, E∞
        push!(get!(levels, i[6], NTuple{7, Float64}[]),
              (i[5], r[1], r[2], r[4], i[1], i[2], i[3]))
    end
    foreach(v -> sort!(v, by = first), values(levels))

    # atomic masses: element records are [n_ions, Z], reals [abundance, mass]
    mass = Dict{Int, Float64}()
    for j in 1:N
        ptr[2, j] == 13 && (mass[ints(j)[2]] = reals(j)[2])
    end

    cand = [j for j in 1:N if ptr[2, j] == type]
    type in RADIATION_TYPES && (cand = [j for j in cand if ptr[5, j] <= MAX_TABLE])   # keep the fixtures small
    rng = MersenneTwister(type)
    if !isempty(selection)
        want = [parse.(Int, split(s, ":")) for s in split(selection, ",")]
        sel = [j for j in cand if (ints(j)[end], ints(j)[7], ints(j)[5]) in Tuple.(want)]
    elseif type in PARENT_TYPES
        # half of the records leave the ground level of the next ion, half an excited one
        kp = parent_positions(type)[1]
        ground = [j for j in cand if ints(j)[kp] == 1]
        excited = [j for j in cand if ints(j)[kp] > 1]
        pick(v, n) = v[randperm(rng, length(v))[1:min(n, length(v))]]
        sel = vcat(pick(ground, nrec ÷ 2), pick(excited, nrec - nrec ÷ 2))
    else
        sel = cand[randperm(rng, length(cand))[1:min(nrec, length(cand))]]
    end
    mkpath(outdir)
    open(joinpath(outdir, "type$(lpad(type, 2, '0'))$tag.in"), "w") do io
        for j in sel
            iv, rv = ints(j), reals(j)
            lv = get(levels, iv[end], NTuple{7, Float64}[])
            for (ci, (t, xpx, xee, f0, f1)) in enumerate(CONDITIONS)
                println(io, type, " ", ptr[3, j], " ", length(rv), " ", length(iv))
                println(io, join(iv, " "))
                println(io, join(rv, " "))
                amass = type == 50 ? get(mass, iv[3], 1.0) : 1.0       # type 50 ints: i, k, Z, ionN
                cfrac = type == 50 ? CFRAC[mod1(ci, length(CFRAC))] : 0.0
                println(io, join((t, xpx, xee, xpx*f0, xpx*f1, 1.0, -1.0, cfrac, 1e-3, 1e-3, 0.5, 0.5, amass), " "))
                # radiation (only the photoionization types need one) and the level of the
                # next ion that a photoionization record leaves
                kpar, gpar, epar = 1, 1.0, 0.0
                if type in PARENT_TYPES
                    ik, ip = parent_positions(type)
                    kpar = iv[ik]                     # parent level of ParPhotoIonize1/2/3
                    pl = get(levels, iv[ip], NTuple{7, Float64}[])
                    if 1 <= kpar <= length(pl)
                        gpar, epar = pl[kpar][3], pl[kpar][2]
                    end
                end
                norm = type in RADIATION_TYPES ? 1e16 : 0.0
                lfast = type in RADIATION_TYPES ? LFAST[mod1(ci, length(LFAST))] : 1
                println(io, join((0.1, 0.0015, 1.0, norm, kpar, gpar, epar, lfast), " "))
                # nlev, then only the levels the record refers to
                need = sort(unique(filter(i -> 1 <= i <= length(lv), vcat(iv[1:end-1], 1, length(lv)))))   # incl. the ground (1) and the continuum (nlev) levels
                println(io, length(lv), " ", length(need))
                for i in need
                    (_, E, g, Einf, n, s2, L) = lv[i]
                    println(io, join((i, E, g, Einf, Int(n), Int(s2), Int(L)), " "))
                end
            end
        end
    end
    println("type $type: ", length(sel), " records x ", length(CONDITIONS), " conditions")
end

main(ARGS[1], ARGS[2], parse(Int, ARGS[3]), length(ARGS) > 3 ? parse(Int, ARGS[4]) : 25,
     length(ARGS) > 4 ? ARGS[5] : "", length(ARGS) > 5 ? ARGS[6] : "")

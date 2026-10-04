# The ionization balance of a mixture of elements: the level balance of each element (`Elements`) at the electron
# density that all of them together produce. XSTAR's calc_hmc_all (the electrons of each element) and dsec (the iteration of
# the electron fraction).

const minimum_abundance = 1e-24         # elements with a smaller abundance (relative to hydrogen) are left out (calc_hmc_all)
const electron_factor = 1.2             # dsec: step of the electron fraction while the solution is not bracketed
const electron_tolerance = 1e-4         # dsec's crite: |xee - electrons|/xee at convergence
const electron_iterations = 100         # at most this many evaluations of the elements in the iteration
const hydrogen_Z = 1

"""
    Mixture(records, levels; multiplier=Dict())

The elements of the database `records` (an `Elements` for each `Atom` record) with their abundances relative to
hydrogen, the `abundance` of the record times `multiplier[Z]` (1 if missing; XSTAR's `habund`, `heabund`, ...). Elements
whose abundance is below `minimum_abundance` (0 for instance) are left out, as in XSTAR. `levels` is `levels(records)`.
"""
struct Mixture{L<:Levels}
    levels::L
    Z::Vector{Int}                  # atomic numbers of the elements
    abundance::Vector{Float64}      # their abundances relative to hydrogen
    elements::Vector{Elements{L}}
end

function Mixture(records, levels::Levels; multiplier=Dict{Int, Float64}())
    atoms = sort!(filter(r -> r isa Atom, records); by=atom -> atom.Z)
    atoms = [atom for atom in atoms if Float64(atom.abundance)*get(multiplier, Int(atom.Z), 1.0) > minimum_abundance]
    Z = [Int(atom.Z) for atom in atoms]
    abundance = [Float64(atom.abundance)*get(multiplier, Int(atom.Z), 1.0) for atom in atoms]
    Mixture(levels, Z, abundance, [Elements(records, levels, z) for z in Z])
end

"""
    electrons(mixture, fractions)

The number of electrons per hydrogen nucleus of the gas with the ion `fractions` of each element (each from the
neutral up and the bare nucleus last, see `ion_fractions`): XSTAR's `enelec`. Every ion of charge `q` has `q`
electrons free and the bare nucleus `Z`; the abundance weights the elements.
"""
function electrons(mixture::Mixture, fractions)
    sum(enumerate(mixture.Z)) do (k, Z)
        f = fractions[k]
        bound = sum(f[1:end - 1])
        mixture.abundance[k]*(sum(f[i]*(i - 1) for i in 1:length(f) - 1) + max(0.0, 1 - bound)*Z)
    end
end

"""
    ionization_balance(mixture, T, ntot; radiation=NO_RADIATION, escape=nothing, xee=1.0, iterate=true,
                       neutral=0.0, lfast=photoionization_lfast, tolerance=electron_tolerance)

The level populations of all the elements of `mixture` at the temperature `T` (10⁴ K) and the hydrogen density
`ntot` (cm⁻³), in the radiation field `radiation`. The electron density is `ntot xee`, and `xee` is iterated until it is
the number of electrons per hydrogen nucleus that the balance gives, to the relative `tolerance`, in XSTAR's way
(`dsec`): the solution is bracketed by steps of a factor 1.2 from the starting value `xee`, then found by false
position. With `iterate=false` the balance is made at `xee` (the reference run of XSTAR does that, with `xee = 1`).
`neutral` is the density of neutral hydrogen (cm⁻³) that the charge exchange rates see at the first evaluation, the
following ones take it from the populations of the hydrogen ground state; `escape` and `lfast` are as for
`element_matrix!`.

Returns a named tuple with the electron fraction `xee` and density `nₑ`, the `electrons` per hydrogen nucleus, the
neutral hydrogen density `nₕ`, the `populations` and the ion `fractions` of each element (in the order of
`mixture.Z`), the number of evaluations `iterations` and whether it `converged`.
"""
function ionization_balance(mixture::Mixture, T, ntot; radiation=NO_RADIATION, escape=nothing, xee=1.0, iterate=true,
        neutral=0.0, lfast=photoionization_lfast, tolerance=electron_tolerance)
    elements = mixture.elements
    matrices = [el.N >= sparse_from ? sparse_matrix(el) : zeros(el.N, el.N) for el in elements]
    populations = [zeros(el.N) for el in elements]
    fractions = [zeros(length(el.ions) + 1) for el in elements]
    nₕ = Float64(neutral)
    hydrogen = findfirst(==(hydrogen_Z), mixture.Z)

    # the balance of every element at the electron fraction x; returns XSTAR's elcter = x - electrons
    function evaluate(x)
        cell = Cell(Float64(T), nₕ, ntot*x, Float64(ntot))
        Threads.@threads for k in eachindex(elements)
            element_matrix!(matrices[k], elements[k], cell; radiation, escape, lfast)
            populations[k] = level_populations(matrices[k])
            fractions[k] = ion_fractions(populations[k], elements[k])
        end
        hydrogen === nothing || (nₕ = ntot*mixture.abundance[hydrogen]*populations[hydrogen][1])
        x - electrons(mixture, fractions)
    end

    iterations = 1
    xee = Float64(xee)
    elcter = evaluate(xee)
    converged = !iterate || abs(elcter)/max(1e-48, xee) < tolerance
    if iterate && !converged
        xee, elcter, iterations, converged = bracket_electron_fraction(evaluate, xee, elcter; tolerance)
    end
    (; xee, nₑ=ntot*xee, electrons=xee - elcter, nₕ, populations, fractions, iterations, converged)
end

# dsec's iteration of the electron fraction: `evaluate(x)` is x minus the electrons that the balance at x gives. A negative
# value means that x is too small. The solution is bracketed by steps of `electron_factor` and then found by false position,
# the evaluation at the returned fraction being the last one made (so that the populations are those of the solution).
function bracket_electron_fraction(evaluate, xee, elcter; tolerance)
    xeel, xeeh = 0.0, 1.0                # the fractions that bound the solution, and the values of elcter there
    elctrl, elctrh = 1.0, -1.0
    ilx = ihx = false
    iterations = 1
    while iterations < electron_iterations
        if elcter < 0
            ihx, xeeh, elctrh = true, xee, elcter
            if !ilx
                xee *= electron_factor
                elcter = evaluate(xee); iterations += 1
                abs(elcter)/max(1e-48, xee) < tolerance && return (xee, elcter, iterations, true)
                continue
            end
        else
            ilx, xeel, elctrl = true, xee, elcter
            if !ihx
                xee /= electron_factor
                elcter = evaluate(xee); iterations += 1
                abs(elcter)/max(1e-48, xee) < tolerance && return (xee, elcter, iterations, true)
                continue
            end
        end
        xee = (xeel*elctrh - xeeh*elctrl)/(elctrh - elctrl)
        elcter = evaluate(xee); iterations += 1
        abs(elcter)/max(1e-48, xee) < tolerance && return (xee, elcter, iterations, true)
    end
    (xee, elcter, iterations, false)
end

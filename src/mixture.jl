# The ionization balance of a mixture of elements: the level balance of each element (`Elements`) at the electron
# density that all of them together produce. XSTAR's calc_hmc_all (the electrons of each element) and dsec (the iteration of
# the electron fraction).

const minimum_abundance = 1e-24         # elements with a smaller abundance (relative to hydrogen) are left out (calc_hmc_all)
const electron_factor = 1.2             # dsec: step of the electron fraction while the solution is not bracketed
const electron_tolerance = 1e-4         # dsec's crite: |xee - electrons|/xee at convergence
const electron_iterations = 100         # at most this many evaluations of the elements in the iteration
const hydrogen_Z = 1
const temperature_factor = 1.2          # dsec: step of the temperature while the solution is not bracketed
const imbalance_tolerance = 1e-4        # dsec's crith: |hmctot| at convergence (the relative heating - cooling)
const temperature_tolerance = 2e-9      # dsec's critt: the iteration stops if the temperature changes by less than this
const large_imbalance = 0.9             # dsec: a second step of the temperature if |hmctot| is larger than this
const temperature_iterations = 99       # at most this many evaluations of the temperature (XSTAR's niter)
const default_T_min = 1e-3              # the lowest temperature of the iteration (10⁴ K)

"""
    gas_density(ntot, T)

The hydrogen density (cm⁻³) of the gas at the temperature `T` (10⁴ K): `ntot` if it is a number (constant density), or `ntot(T)` if it is a function, such as `ConstantPressure(P)`.
"""
gas_density(ntot::Number, T) = ntot
gas_density(ntot, T) = ntot(T)

"""
    gas_density(ntot, T, r)

The same at the distance `r` (cm) from the source: only `PowerLawDensity` depends on it.
"""
gas_density(ntot, T, r) = gas_density(ntot, T)

const pressure_per_density = Float64(1.38f-12)    # `xpx = p/1.38e-12/t`: k × 10⁴ K, as a single-precision literal
const pressure_T_floor = 1e-24

"""
    ConstantPressure(P)

The hydrogen density at the constant pressure `P` (dyn cm⁻²), as a function of the temperature `T` (10⁴ K): `P/(1.38×10⁻¹² T)` (XSTAR's `lcpres=1`, in which the pressure is that of the hydrogen, `n_H k T`).
"""
struct ConstantPressure
    P::Float64
end
(p::ConstantPressure)(T) = p.P/pressure_per_density/max(T, pressure_T_floor)

"""
    PowerLawDensity(n0, r0, exponent)

The hydrogen density `n0 (r/r0)^exponent` at the distance `r` of the source (XSTAR's `radexp`, for the constant density option): `gas_density(d, T, r)`; without `r` it is `n0`.
"""
struct PowerLawDensity
    n0::Float64
    r0::Float64
    exponent::Float64
end
gas_density(d::PowerLawDensity, T) = d.n0
gas_density(d::PowerLawDensity, T, r) = d.n0*(r/d.r0)^d.exponent

"""
    Mixture(records, levels; abundances=nothing, multiplier=Dict())

The elements of the database `records` (an `Elements` for each `Atom` record) with their abundances relative to
hydrogen, the `abundance` of the record (or `abundances[Z]`, a table of the elements 1 to 30 such as `abundance_table(:angr)`) times `multiplier[Z]` (1 if missing; XSTAR's `habund`, `heabund`, ...). Elements
whose abundance is below `minimum_abundance` (0 for instance) are left out, as in XSTAR. `levels` is `levels(records)`.
"""
struct Mixture{L<:Levels}
    levels::L
    Z::Vector{Int}                  # atomic numbers of the elements
    abundance::Vector{Float64}      # their abundances relative to hydrogen
    elements::Vector{Elements{L}}
end

function Mixture(records, levels::Levels; abundances=nothing, multiplier=Dict{Int, Float64}())
    atoms = sort!(filter(r -> r isa Atom, records); by=atom -> atom.Z)
    relative(atom) = (abundances === nothing ? Float64(atom.abundance) : Float64(abundances[Int(atom.Z)]))*get(multiplier, Int(atom.Z), 1.0)
    atoms = [atom for atom in atoms if relative(atom) > minimum_abundance]
    Z = [Int(atom.Z) for atom in atoms]
    abundance = [relative(atom) for atom in atoms]
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

# `escape` of one element: a vector with one vector of pairs for each element (see `escape_probabilities`) gives that of the element k;
# anything else (nothing, a function of the record, a vector of pairs of the records) is the same for all
element_escape(escape, k) = escape isa AbstractVector && !isempty(escape) && first(escape) isa AbstractVector ? escape[k] : escape

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
    memos = [RatesMemo(el) for el in elements]       # (the cross sections are integrated once for all the electron fractions)
    nₕ = Float64(neutral)
    hydrogen = findfirst(==(hydrogen_Z), mixture.Z)

    # the balance of every element at the electron fraction x; returns XSTAR's elcter = x - electrons
    function evaluate(x)
        cell = Cell(Float64(T), nₕ, ntot*x, Float64(ntot))
        Threads.@threads for k in eachindex(elements)
            element_matrix!(matrices[k], elements[k], cell; radiation, escape=element_escape(escape, k), lfast, memo=memos[k])
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

function process_term(process, radiation, T, nₑ)
    h, c = heating_cooling(process, radiation, T, nₑ)
    (; process, heating=h, cooling=c)
end

"""
    heating_cooling(mixture, balance, T, ntot, processes; radiation=NO_RADIATION, escape=nothing, lfast=photoionization_lfast)

The heating and cooling (erg cm⁻³ s⁻¹) of the gas of `mixture` at the temperature `T` (10⁴ K) and hydrogen density `ntot`
whose populations are those of `balance` (from `ionization_balance`, at the same `T`, `ntot` and `radiation`). `processes` is a collection of
`AbstractContinuum` processes (`standard_processes(compton)`: Compton and Thomson scattering, free-free absorption and bremsstrahlung): the totals of XSTAR's `calc_hmc_all` and
`heatf` are the elements (`element_heating`, weighted by their abundances) plus the `heating` and `cooling` derived from each process.

Returns a named tuple with `heating` and `cooling` (the radiative energy that the gas absorbs and emits: `httot`, `cltot`), `imbalance` =
`2 (heating - cooling)/(heating + cooling)` (`hmctot`, what the temperature is iterated to zero), `heating2` and `cooling2` (the same for
the energy of the electrons), `processes`, a vector with the `process`, its `heating` and its `cooling`, and `elements`, the `heating`, `cooling`,
`heating2` and `cooling2` of each element weighted by its abundance, in the order of `mixture.Z`.
"""
function heating_cooling(mixture::Mixture, balance, T, ntot, processes; radiation=NO_RADIATION, escape=nothing,
        lfast=photoionization_lfast)
    cell = Cell(Float64(T), Float64(balance.nₕ), Float64(balance.nₑ), Float64(ntot))
    elements = Vector{NTuple{4, Float64}}(undef, length(mixture.elements))
    Threads.@threads for k in eachindex(mixture.elements)
        elements[k] = element_heating(mixture.elements[k], cell, balance.populations[k]; radiation, escape=element_escape(escape, k), lfast)
    end
    weighted = [mixture.abundance[k] .* elements[k] for k in eachindex(elements)]
    terms = [process_term(process, radiation, T, balance.nₑ) for process in processes]
    process_heating, process_cooling = sum(t.heating for t in terms; init=0.0), sum(t.cooling for t in terms; init=0.0)
    total_heating = sum(w[1] for w in weighted) + process_heating
    total_cooling = sum(w[2] for w in weighted) + process_cooling
    (; heating=total_heating, cooling=total_cooling, imbalance=2*(total_heating - total_cooling)/(heating_floor + total_heating + total_cooling),
       heating2=sum(w[3] for w in weighted) + process_heating, cooling2=sum(w[4] for w in weighted) + process_cooling,
       processes=terms,
       elements=[(; heating=w[1], cooling=w[2], heating2=w[3], cooling2=w[4]) for w in weighted])
end

"""
    thermal_equilibrium(mixture, ntot, processes; radiation=NO_RADIATION, escape=nothing, T=1.0, xee=1.0, T_min=default_T_min,
                        tolerance=imbalance_tolerance, iterations=temperature_iterations)

The temperature at which the heating and cooling of the gas balance, and the ionization balance there, for the hydrogen
density `ntot` and the radiation field `radiation`, starting from the temperature `T` (10⁴ K) and the electron fraction `xee`
and iterating as XSTAR's `dsec` does: at each temperature the electron fraction is iterated (`ionization_balance`) and the
heating and cooling summed (`heating_cooling`); the temperature is changed by factors of 1.2 (of 1.44 when the relative
imbalance exceeds 0.9) until it brackets the solution, then found by false position, with the stale end halved, until the
relative imbalance `|2 (heating - cooling)/(heating + cooling)|` is below `tolerance` (or the temperature stops changing,
or `iterations` temperatures were tried). `T_min` (10⁴ K) is the lowest temperature.

Returns the named tuple of `heating_cooling` for the last temperature with the added fields `T`, `xee`, `nₑ`, `nₕ`,
`populations`, `fractions` (those of `ionization_balance`), `evaluations` (the temperatures tried), `converged`,
and `trace`, the temperatures tried with their imbalance.
"""
function thermal_equilibrium(mixture::Mixture, ntot, processes; radiation=NO_RADIATION, escape=nothing, T=1.0,
        xee=1.0, T_min=default_T_min, lfast=photoionization_lfast, tolerance=imbalance_tolerance, iterations=temperature_iterations)
    xee = Float64(xee)
    neutral = 0.0
    local balance, energy
    function evaluate(t)
        n = gas_density(ntot, t)
        balance = ionization_balance(mixture, t, n; radiation, escape, xee, neutral, lfast)
        energy = heating_cooling(mixture, balance, t, n, processes; radiation, escape, lfast)
        xee, neutral = balance.xee, balance.nₕ
        energy.imbalance
    end
    t, converged, evaluations, trace = bracket_temperature(evaluate, Float64(T); tolerance, iterations, T_min)
    (; energy..., T=t, ntot=gas_density(ntot, t), xee=balance.xee, nₑ=balance.nₑ, nₕ=balance.nₕ, populations=balance.populations, fractions=balance.fractions,
       evaluations, converged, trace)
end

# dsec's iteration of the temperature: `evaluate(t)` is the relative imbalance of the heating and cooling at the temperature t,
# negative when the cooling exceeds the heating (t too high). The solution is bracketed by steps of `temperature_factor` (twice
# that if the imbalance exceeds `large_imbalance`) and then found by false position, with the end of the bracket that has not
# changed halved, as XSTAR does. Returns the last temperature, whether `|imbalance| <= tolerance` there, the number of
# evaluations and the trace of (temperature, imbalance).
function bracket_temperature(evaluate, t; tolerance=imbalance_tolerance, iterations=temperature_iterations, T_min=default_T_min)
    t = max(t, T_min)
    tl = th = hmcttl = hmctth = 0.0
    to = 1e30
    iht = ilt = iuht = iult = false
    trace = Tuple{Float64, Float64}[]
    nnt = 0
    while true
        hmctot = evaluate(t)
        nnt += 1
        push!(trace, (t, hmctot))
        abs(hmctot) <= tolerance && break
        nnt < iterations || break
        if hmctot < 0                                       # the cooling exceeds the heating: the temperature is too high
            iht, th, hmctth, iuht = true, t, hmctot, true
            iult || (hmcttl /= 2)
            iult = false
            if !ilt
                t /= temperature_factor
                abs(hmctot) > large_imbalance && (t /= temperature_factor)
                t = max(t, T_min)
                continue
            end
        else
            ilt, tl, hmcttl, iult = true, t, hmctot, true
            iuht || (hmctth /= 2)
            iuht = false
            if !iht
                t *= temperature_factor
                abs(hmctot) > large_imbalance && (t *= temperature_factor)
                continue
            end
        end
        abs(1 - t/to) < temperature_tolerance && break      # (not converging)
        to = t
        t = (tl*hmctth - th*hmcttl)/(hmctth - hmcttl)
    end
    (t, abs(trace[end][2]) <= tolerance, nnt, trace)
end

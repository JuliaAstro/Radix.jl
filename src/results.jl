# The results of a run as arrays, and their single file: the zone table, the spectrum and the tables of lines and edges (XSTAR's xout_abund1, xout_cont1/xout_spect1, xout_lines1 and xout_rrc1).

const roman_numerals = ["i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x", "xi", "xii", "xiii", "xiv", "xv", "xvi", "xvii", "xviii", "xix", "xx",
    "xxi", "xxii", "xxiii", "xxiv", "xxv", "xxvi", "xxvii", "xxviii", "xxix", "xxx", "xxxi"]
const results_min_fraction = 1e-15          # binemis: the lines with less than this fraction of the luminosity of the source are not in the tables
const band_E_min = 13.6                     # the band of the luminosity of the source (eV)
const band_E_max = 1.36e4

# `he_ii` → `he_iii`: the label of the ion that a bare nucleus is
function next_stage(label)
    element, stage = rsplit(label, '_'; limit=2)
    n = findfirst(==(stage), roman_numerals)
    n === nothing ? label * "+" : element * "_" * roman_numerals[n + 1]
end

# the threshold (eV) of a recombination edge: as `rank_energy` has it, else the ionization energy of its level (the cross section of some records starts at 0 above the threshold), else NaN
function threshold_energy(c)
    E = rank_energy(c)
    E > 0 && return E
    level = hasproperty(c, :level) ? get(c.levels, (c.ion, c.level), nothing) : nothing
    level === nothing && return NaN
    Float64(level.E_inf) - Float64(level.E)
end

# the luminosity (erg s⁻¹) of the incident spectrum in the band of the source
band_luminosity(E, incident) = sum((incident[i] + incident[i - 1])*(E[i] - E[i - 1])/2 for i in 2:length(E) if band_E_min <= E[i] <= band_E_max; init=0.0)*constants().ergsev

"""
    Results

The tables of a run (`tables`): the fields `zones`, `spectrum`, `lines` and `edges` are named tuples of columns. `write(path, results)` writes them to one FITS file.
"""
struct Results{Z, S, L, E}
    zones::Z
    spectrum::S
    lines::L
    edges::E
end

"""
    tables(model, run; min_fraction=1e-15)

The results of `slab_model(model; ...)` as `Results`, four tables (their fields), each a named tuple of columns (vectors; units below):

- `zones`: one row for each zone (the state at its inner edge): `zone`, `radius` (cm), `depth` (cm, of the inner edge from the first), `thickness` (cm), `ion_parameter` (`L/(n r²)`), `x_e`, `n_h` and `n_e` (cm⁻³),
  `pressure` (dyn cm⁻²), `temperature` (10⁴ K), `heating` and `cooling` (erg cm⁻³ s⁻¹), `heat_error` (relative) and the fraction of every ion, named by its label (`h_i`, `he_ii`, ..., with the bare nuclei);
- `spectrum`: for each energy, `energy` (eV), `incident`, `transmitted` (with the continuum opacity: the transmitted column of XSTAR), `transmitted_lines` (with the lines and free-free too, as the spectrum is attenuated along the slab) and `emit_inward` and `emit_outward`
  (the diffuse emission of the continuum, summed over the zones), all in erg s⁻¹ erg⁻¹;
- `lines`: the lines (rate type 4) with more than `min_fraction` of the luminosity of the source or a depth: `ion`, `lower`, `upper` (labels of the levels), `wavelength` (Å), `emit_inward` and `emit_outward` (erg s⁻¹), `depth_inward` and `depth_outward` (at the centre of the line);
- `edges`: the same for the recombination continua (rate type 7), with the threshold `energy` (eV) in place of the wavelength.

`inward` is the first direction of the escape probabilities (XSTAR's `emit_inward`). Arrays only: `write(path, results)` writes them to one file.
"""
function tables(model::Model, run; min_fraction=results_min_fraction)
    mixture, labels = model.mixture, model.ion_labels
    K = constants()
    zones = run.zones
    luminosity = band_luminosity(run.E, run.incident)
    depth = [0.0; cumsum([z.Δr for z in zones])[1:end - 1]]
    fractions = Pair{Symbol, Vector{Float64}}[]
    for (k, element) in enumerate(mixture.elements)
        names = [labels[Int(ion)] for ion in element.ions]
        push!(names, next_stage(names[end]))
        for (i, name) in enumerate(names)
            push!(fractions, Symbol(name) => [z.fractions[k][i] for z in zones])
        end
    end
    zone_table = (; zone=collect(1:length(zones)), radius=[z.r for z in zones], depth, thickness=[z.Δr for z in zones],
        ion_parameter=[luminosity/(z.ntot*z.r^2) for z in zones], x_e=[z.xee for z in zones], n_h=[z.ntot for z in zones], n_e=[z.nₑ for z in zones],
        pressure=[z.ntot*pressure_per_density*z.T for z in zones], temperature=[z.T for z in zones], heating=[z.heating for z in zones],
        cooling=[z.cooling for z in zones], heat_error=[z.imbalance for z in zones], fractions...)
    spectrum = (; energy=run.E, incident=run.incident, transmitted=run.transmitted, transmitted_lines=run.incident .* exp.(-run.dpthc),
        emit_inward=run.spectrum.inward_continuum, emit_outward=run.spectrum.outward_continuum)
    levels = mixture.levels
    label(ion, level) = String(strip(levels[(Int(ion), Int(level))].label))
    threshold = min_fraction*luminosity
    line_rows = NamedTuple[]
    edge_rows = NamedTuple[]
    for (k, element) in enumerate(mixture.elements), j in eachindex(element.rates)
        c = element.rates[j]
        inward, outward = run.luminosities.inward[k][j]*luminosity_unit, run.luminosities.outward[k][j]*luminosity_unit
        if c isa AtomicLine2 && c.rtype == line_data_type && element.lo[j] != 0 && inward + outward > threshold
            push!(line_rows, (; ion=labels[Int(c.ion)], lower=label(c.ion, c.transition.lower), upper=label(c.ion, c.transition.upper), wavelength=abs(Float64(c.λ)),
                emit_inward=inward, emit_outward=outward, depth_inward=run.depths.inward[k][j], depth_outward=run.depths.outward[k][j]))
        elseif c isa opacity_edge_types && c.rtype == edge_rate_type && element.lo[j] != 0 && inward + outward > threshold
            push!(edge_rows, (; ion=labels[Int(c.ion)], level=label(c.ion, c.level), energy=threshold_energy(c), emit_inward=inward, emit_outward=outward,
                depth_inward=run.depths.inward[k][j], depth_outward=run.depths.outward[k][j]))
        end
    end
    columns(rows, names) = NamedTuple{names}(Tuple([row[name] for row in rows] for name in names))
    Results(zone_table, spectrum,
        columns(line_rows, (:ion, :lower, :upper, :wavelength, :emit_inward, :emit_outward, :depth_inward, :depth_outward)),
        columns(edge_rows, (:ion, :level, :energy, :emit_inward, :emit_outward, :depth_inward, :depth_outward)))
end

const column_units = Dict(:radius => "cm", :depth => "cm", :thickness => "cm", :n_h => "cm**-3", :n_e => "cm**-3", :pressure => "dyn/cm**2", :temperature => "1e4 K",
    :heating => "erg/cm**3/s", :cooling => "erg/cm**3/s", :energy => "eV", :incident => "erg/s/erg", :transmitted => "erg/s/erg", :transmitted_lines => "erg/s/erg",
    :emit_inward => "erg/s", :emit_outward => "erg/s", :wavelength => "Angstrom")
const spectrum_units = Dict(:emit_inward => "erg/s/erg", :emit_outward => "erg/s/erg")      # (the diffuse emission of the spectrum is per unit of energy)

# a column that FITS can hold: numbers as they are, labels padded to the same length
fits_column(v::AbstractVector{<:AbstractString}) = (width = maximum(length, v; init=1); [rpad(s, width) for s in v])
fits_column(v::AbstractVector{<:Integer}) = Int64.(v)
fits_column(v::AbstractVector{<:AbstractFloat}) = Float64.(v)

"""
    write(path, results::Results; meta=(;), overwrite=true)

Writes the `results` of `tables(model, run)` to one FITS file: an empty primary HDU whose header has the `meta` (a named tuple of keyword values, e.g. the parameters of the run) and a binary table extension for
each of `ZONES`, `SPECTRUM`, `LINES` and `EDGES`, with the units of the columns.
"""
function Base.write(path::AbstractString, results::Results; meta=(;), overwrite=true)
    !overwrite && isfile(path) && throw(ArgumentError("$path exists"))
    deck() = Vector{FITSFiles.Card{<:Any}}()
    cards = deck()
    push!(cards, FITSFiles.Card("CREATOR", "Radix.jl"))
    for (key, value) in pairs(meta)
        push!(cards, FITSFiles.Card(uppercase(String(key)), value))
    end
    hdus = FITSFiles.HDU[FITSFiles.HDU(FITSFiles.Primary, missing, cards)]
    for (name, table) in zip((:zones, :spectrum, :lines, :edges), (results.zones, results.spectrum, results.lines, results.edges))
        data = NamedTuple{keys(table)}(Tuple(fits_column(v) for v in values(table)))
        cards = deck()
        push!(cards, FITSFiles.Card("EXTNAME", uppercase(String(name))))
        for (i, column) in enumerate(keys(table))
            unit = name === :spectrum ? get(spectrum_units, column, get(column_units, column, nothing)) : get(column_units, column, nothing)
            unit === nothing || push!(cards, FITSFiles.Card("TUNIT$i", unit))
        end
        push!(hdus, FITSFiles.HDU(FITSFiles.Bintable, data, cards))
    end
    Base.write(path, hdus)
    path
end

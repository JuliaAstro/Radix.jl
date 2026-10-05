# Heating and cooling of the gas: the energies of the level populations, Compton, free-free and bremsstrahlung, and the
# temperature at which they balance.

using FITSFiles

# the term of a process in the result of heating_cooling
term(hc, P) = only(t for t in hc.processes if t.process isa P)

function toy_heating_tests()
    @testset "Heating and cooling" begin
        K = Radix.constants()

        @testset "the Compton table" begin
            # coheat.dat: i j sx x e d (the third, fifth and sixth columns are read)
            text = join(["    $i   $j  $(sx)E+00  0.1000E+08  $(e)E+00  $(round(1.0 + 2.0*sx + 3.0*e, digits=4))E+00\n"
                         for (i, sx) in enumerate((1.0, 2.0, 4.0)) for (j, e) in enumerate((0.5, 1.0, 3.0))], "")
            table = Radix.load(Radix.Compton, IOBuffer(text))
            @test table.Te == [1.0, 2.0, 4.0] && table.Eph == [0.5, 1.0, 3.0] && size(table.Es) == (3, 3)
            @test table.Es[2, 3] == 1.0 + 2.0*2.0 + 3.0*3.0
            # from a file, as the stream
            mktemp() do path, io
                write(io, text); close(io)
                from_file = Radix.load(Radix.Compton, path)
                @test from_file.Te == table.Te && from_file.Eph == table.Eph && from_file.Es == table.Es
            end
            # interpolation of a linear function is exact, below the limit energy the limit 4 Te - Eph is used
            for (Eph, Te) in ((0.7, 1.5), (2.0, 3.0), (1.2, 2.2))
                @test Radix.σ(table, Eph, Te) ≈ 1.0 + 2.0*Te + 3.0*Eph
            end
            @test Radix.σ(table, 5e-5, 0.3) == 4*0.3 - 5e-5
        end

        @testset "the integrals" begin
            table = Radix.load(Radix.Compton, IOBuffer(join(["    $i   $j  $(sx)E+00  0.1000E+08  $(e)E+00  $(1.0 + 2.0*sx + 3.0*e)E+00\n"
                for (i, sx) in enumerate((1e-3, 1e-2, 1e-1)) for (j, e) in enumerate((1e-3, 1e-1, 10.0))], "")))
            E = [1e3, 1e4, 1e5]; F = [3.0, 2.0, 1.0]
            T = 1.0
            ekt = T*K.kT_eV
            Te = (ekt + Radix.compton_Te_floor)/K.electron_rest_eV
            cmp1, cmp2 = Radix.integral(table, Radix.Radiation(E, F), T)
            # by hand: heating σ ∫ F E/m c² dE, and the cooling coefficient (σ ∫ F cmpfnc dE + heating)/kT
            trap(g) = sum((g[k] + g[k - 1])*(E[k] - E[k - 1])/2 for k in 2:3)
            hfake = K.sigma_thomson*trap(F .* E ./ K.electron_rest_eV)
            cohc = -K.sigma_thomson*trap(F .* Radix.σ.(Ref(table), E ./ K.electron_rest_eV, Te))
            @test cmp1 ≈ hfake
            @test cmp2 ≈ (-cohc + hfake)/ekt
            # the process: heating cmp1 nₑ, cooling kT cmp2 nₑ (eV → erg), from one pass for both
            compton = table
            @test compton isa Radix.AbstractContinuum
            @test Radix.heating(compton, Radix.Radiation(E, F), T, 1e4) ≈ cmp1*1e4*K.ergsev
            @test Radix.cooling(compton, Radix.Radiation(E, F), T, 1e4) ≈ ekt*cmp2*1e4*K.ergsev
            @test Radix.heating_cooling(compton, Radix.Radiation(E, F), T, 1e4) ==
                  (Radix.heating(compton, Radix.Radiation(E, F), T, 1e4), Radix.cooling(compton, Radix.Radiation(E, F), T, 1e4))
        end

        @testset "free-free and bremsstrahlung" begin
            E = [1e2, 1e3, 1e4]; F = [5.0, 4.0, 3.0]
            rad = Radix.Radiation(E, F)
            T, nₑ = 2.0, 1e4
            ekt = T*K.kT_eV
            α(e) = K.ff_absorption_coeff*nₑ*(Radix.ion_density_factor*nₑ)/sqrt(T)/e^3*(1 - exp(-e/ekt))
            j(e) = K.ff_emission_coeff*nₑ*(Radix.ion_density_factor*nₑ)*exp(-e/ekt)/sqrt(T)
            @test Radix.heating(Radix.FreeFree(), rad, T, nₑ) ≈ K.ergsev*sum((F[k]*α(E[k]) + F[k - 1]*α(E[k - 1]))*(E[k] - E[k - 1])/2 for k in 2:3)
            @test Radix.cooling(Radix.Bremsstrahlung(), rad, T, nₑ) ≈ K.ergsev*sum((j(E[k]) + j(E[k - 1]))*(E[k] - E[k - 1])/2 for k in 2:3)
            @test Radix.cooling(Radix.Bremsstrahlung(), rad, T, 2nₑ) ≈ 4*Radix.cooling(Radix.Bremsstrahlung(), rad, T, nₑ)       # ∝ nₑ²
            @test Radix.heating(Radix.FreeFree(), Radix.Radiation(E, 2 .* F), T, nₑ) ≈ 2*Radix.heating(Radix.FreeFree(), rad, T, nₑ)
            # a process that does not heat or cool has 0; heating_cooling gives both
            @test Radix.cooling(Radix.FreeFree(), rad, T, nₑ) == 0 && Radix.heating(Radix.Bremsstrahlung(), rad, T, nₑ) == 0
            @test Radix.heating_cooling(Radix.FreeFree(), rad, T, nₑ) == (Radix.heating(Radix.FreeFree(), rad, T, nₑ), 0.0)
            @test Radix.heating_cooling(Radix.Bremsstrahlung(), rad, T, nₑ) == (0.0, Radix.cooling(Radix.Bremsstrahlung(), rad, T, nₑ))
            @test all(p -> p isa Radix.AbstractContinuum, (Radix.FreeFree(), Radix.Bremsstrahlung(), Radix.Thomson()))
        end

        @testset "the processes and their opacities" begin
            E = [1e2, 1e3, 1e4]
            rad = Radix.Radiation(E, [5.0, 4.0, 3.0])
            T, nₑ = 2.0, 1e4
            ekt = T*K.kT_eV
            # free-free: the absorption coefficient, with the stimulated emission; in the total opacity only
            ff = Radix.α(Radix.FreeFree(), E, T)
            @test ff ≈ [K.ff_absorption_coeff/sqrt(T)/e^3*(1 - exp(-e/ekt)) for e in E]
            @test Radix.α(Radix.FreeFree(), rad, T) == ff                          # (a Radiation gives its grid)
            # the opacity is the coefficient times the density that absorbs: nₑ nᵢ with nᵢ = 1.4 nₑ
            @test Radix.density(Radix.FreeFree(), nₑ) == nₑ*Radix.ion_density_factor*nₑ
            @test Radix.opacity(Radix.FreeFree(), E, T, nₑ) ≈ Radix.density(Radix.FreeFree(), nₑ) .* ff
            @test !(Radix.FreeFree() isa Radix.AbstractScattering)
            # Thomson scattering: nₑ σ_T (1 - cfrac) at every energy, in the continuum as well
            @test Radix.σ(Radix.Thomson(), E, T) == fill(K.sigma_thomson, 3)               # the coefficient, per electron
            @test Radix.σ(Radix.Thomson(0.25), E, T) ≈ fill(0.75*K.sigma_thomson, 3)
            @test Radix.density(Radix.Thomson(), nₑ) == nₑ
            @test Radix.opacity(Radix.Thomson(0.25), E, T, nₑ) ≈ fill(0.75*nₑ*K.sigma_thomson, 3)
            @test Radix.Thomson().cfrac == 0 && (Radix.Thomson() isa Radix.AbstractScattering)
            # the coefficients from CODATA 2022 agree with those that XSTAR writes (Gaunt factor 1) to 0.1%
            xstar = Radix.ucalc_constants()
            @test K.ff_absorption_coeff ≈ xstar.ff_absorption_coeff rtol=1e-3
            @test K.ff_emission_coeff ≈ xstar.ff_emission_coeff rtol=1e-3
            @test K.ff_absorption_coeff ≈ 2.611798464985841e-37 && K.ff_emission_coeff ≈ 1.0325265054853345e-13
            # the processes without opacity have none, and Thomson scattering neither heats nor cools
            @test Radix.α(Radix.Bremsstrahlung(), E, T) == zeros(3) && Radix.j(Radix.Bremsstrahlung(), E, T) == Radix.j(Radix.Bremsstrahlung(), rad, T) && Radix.opacity(Radix.Bremsstrahlung(), E, T, nₑ) == zeros(3)
            compton0 = Radix.Compton([1e-4, 1.0], [1e-3, 1.0], [1.0 2.0; 3.0 4.0])
            @test Radix.σ(compton0, rad, T) == fill(K.sigma_thomson, 3) && Radix.opacity(compton0, rad, T, nₑ) == zeros(3)
            @test Radix.heating_cooling(Radix.Thomson(), rad, T, nₑ) == (0.0, 0.0)
            # the processes of XSTAR
            table = Radix.Compton([1e-4, 1.0], [1e-3, 1.0], [1.0 2.0; 3.0 4.0])
            processes = Radix.standard_processes(table; cfrac=0.1)
            @test map(typeof, processes) == (Radix.Compton, Radix.FreeFree, Radix.Bremsstrahlung, Radix.Thomson)
            @test processes[1] === table && processes[4].cfrac == 0.1
        end

        @testset "the mapped spectrum" begin
            E = Radix.xstar_energy_grid(999)
            F = 1e10 ./ E
            rad = Radix.Radiation(E, F)
            mapped = Radix.map_spectrum(rad)
            n = 999 - max(2, 999 ÷ 50)                                  # the last bin that nbinc uses
            @test mapped.E === rad.E
            @test mapped.F[1:n - 1] == F[1:n - 1]                       # identical below the top bins
            @test all(mapped.F[n:end] .== F[n])                         # a flat tail above: the flux of the last bin used
            @test mapped.F[end] > F[end]
            # another grid: each bin takes the flux of the nearest bin
            other = Radix.map_spectrum(rad, E[1:10:500])
            @test other.F == F[1:10:500] && length(other.E) == 50
        end

        f32 = Float32
        lv(level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1, 2, 0, 1, level, 5], f32[E, g, 1, 13.6])
        levels = Radix.levels([lv(1, 0.0, 2), lv(2, 10.2, 8), lv(3, 13.6, 1)]; masses=Dict(5 => 1.0))
        cell = Radix.Cell(1.0, 0.0, 1e4, 1e6)
        line = Radix.AtomicLine2(Int32(4), "", Int32[2, 1, 1, 5], f32[1215.67, 0.4162, 6.265e8], levels)
        elements = Radix.Elements([[line]], levels, [5])

        @testset "the energies of a line" begin
            thin = Radix.rate(line, cell)
            # ucalc's ans3 = -(decay energy), ans4 = -(photoexcitation energy)
            @test Radix.ucalc_energies(line, thin) == (-thin.fenergy, -thin.ienergy, 0.0, 0.0)
            # (a type without energies in ucalc has none here)
            @test Radix.ucalc_energies(line, (; fenergy=2.0, ienergy=3.0)) == (-2.0, -3.0, 0.0, 0.0)
        end

        @testset "heating and cooling of the populations" begin
            x = [0.7, 0.3, 0.0]
            ht, cl, ht2, cl2 = Radix.element_heating(elements, cell, x)
            thin = Radix.rate(line, cell)
            # the decay emits energy at the upper level: cooling x₂ × energy of the decay × n; no heating without radiation
            @test ht == 0 && cl ≈ 0.3*thin.fenergy*cell.ntot
            @test ht2 == 0 && cl2 == 0                                  # the line has no energy of the electrons
            # with radiation the photoexcitation absorbs energy at the lower level: a heating x₁ × n × energy of the excitations
            E = Radix.xstar_energy_grid()
            rad = Radix.point_source(E, 1e30 ./ E, 1e13)
            ht, cl, _, _ = Radix.element_heating(elements, cell, x; radiation=rad)
            r = Radix.rate(line, cell; radiation=rad)
            @test r.ienergy > 0
            @test ht ≈ 0.7*r.ienergy*cell.ntot
            @test cl ≈ 0.3*r.fenergy*cell.ntot
            # trapped lines (a smaller escape probability) cool less
            thin_cooling = Radix.element_heating(elements, cell, x)[2]
            trapped_cooling = Radix.element_heating(elements, cell, x; escape=[(0.1, 0.1)])[2]
            @test trapped_cooling ≈ 0.2*thin_cooling
        end

        @testset "the totals" begin
            # no radiation: only the bremsstrahlung (T and nₑ dependent) cools besides the line
            ComptonToy = Radix.Compton([1e-4, 1.0], [1e-3, 1.0], [1.0 2.0; 3.0 4.0])
            mixture = Radix.Mixture(levels, [1], [0.5], [elements])
            balance = (; nₑ=1e4, nₕ=0.0, populations=[[0.7, 0.3, 0.0]])
            processes = Radix.standard_processes(ComptonToy)
            hc = Radix.heating_cooling(mixture, balance, 1.0, 1e6, processes)
            e = Radix.element_heating(elements, Radix.Cell(1.0, 0.0, 1e4, 1e6), [0.7, 0.3, 0.0])
            @test hc.elements[1].cooling ≈ 0.5*e[2] && hc.elements[1].heating == 0
            @test term(hc, Radix.FreeFree).heating == 0 && term(hc, Radix.Compton).heating == 0
            @test length(hc.processes) == 4 && term(hc, Radix.Thomson).heating == 0 && term(hc, Radix.Thomson).cooling == 0
            @test term(hc, Radix.Bremsstrahlung).cooling ≈ Radix.cooling(Radix.Bremsstrahlung(), Radix.NO_RADIATION, 1.0, 1e4)
            @test hc.cooling ≈ 0.5*e[2] + term(hc, Radix.Bremsstrahlung).cooling + term(hc, Radix.Compton).cooling
            @test hc.imbalance ≈ 2*(hc.heating - hc.cooling)/(1e-37 + hc.heating + hc.cooling) && hc.imbalance < 0
        end

        @testset "the temperature iteration" begin
            # the imbalance falls with the temperature, as the cooling grows: found from above and below, the end of the
            # bracket that has not moved is halved in the false position, and a start below the lowest temperature is moved up
            f(t) = (10 - t)/10
            for start in (1.0, 100.0, 10.0, 3.7)
                t, converged, n, trace = Radix.bracket_temperature(f, start)
                @test converged
                @test t ≈ 10 rtol=1e-4
                @test first(trace[1]) == start && n == length(trace) && n <= 30
            end
            # a larger imbalance than 0.9 doubles the step: from T = 1000 the first move is by a factor 1.44
            _, _, _, trace = Radix.bracket_temperature(t -> (10 - t)/10, 1000.0)
            @test trace[2][1] ≈ 1000/1.2^2
            _, _, _, trace = Radix.bracket_temperature(t -> (10 - t)/10, 0.1)
            @test trace[2][1] ≈ 0.1*1.2^2
            # the lowest temperature, and the limit on the number of evaluations
            t, converged, n, _ = Radix.bracket_temperature(t -> -1.0, 1.0; T_min=0.5)
            @test !converged && t == 0.5
            t, converged, n, _ = Radix.bracket_temperature(t -> (10 - t)/10, 1.0; iterations=3)
            @test !converged && n == 3
            # a steeper function (the imbalance of a cooling that grows as T²), and the tolerance
            g(t) = 2*(25 - t^2)/(25 + t^2)
            t, converged, _, _ = Radix.bracket_temperature(g, 2.0; tolerance=1e-8)
            @test converged
            @test t ≈ 5 rtol=1e-6
        end
    end
end

# the reference run of XSTAR (test/reference/xstar_pow_xi2 and xstar_pow_xi2_eq) for all the elements: the heating and cooling of its
# first zone at T = 10⁶ K, and the temperature that it reaches when it iterates
function heating_balance_tests(db)
    @testset "Heating and cooling of the database" begin
        levels = Radix.levels(db)
        coheat = joinpath(dirname(get(ENV, "RADIX_ATDB", "")), "coheat.dat")
        if !isfile(coheat)
            @info "Skipping the heating tests (coheat.dat next to atdb.fits is needed)"
            return
        end
        compton = Radix.load(Radix.Compton, coheat)
        processes = Radix.standard_processes(compton)
        @test size(compton.Es) == (101, 101) && compton.Te[1] == 1e-7 && compton.Eph[end] == 100.0
        mixture = Radix.Mixture(db, levels; multiplier=Dict(3 => 0.0, 4 => 0.0, 5 => 0.0))
        dir = joinpath(@__DIR__, "reference", "xstar_pow_xi2")
        spectrum = fits(joinpath(dir, "xout_cont1.fits"))[3].data
        E, L = Float64.(spectrum.energy), Float64.(spectrum.incident)*1e38
        radius = 1e13                                   # the first zone of the run: log ξ = 2
        incident = Radix.point_source(E, L, radius)
        mapped = Radix.map_spectrum(incident)

        @testset "the Compton integrals of the real comp2" begin
            # drvcomp (test/reference/ucalc) around XSTAR's comp2 on this spectrum at T = 10⁶ K
            exact = Radix.Radiation(E, L ./ (4π*radius^2))             # (the driver's input has F = L/(4π r²), `point_source` the 12.56 of XSTAR)
            cmp1, cmp2 = Radix.integral(compton, exact, 100.0)
            @test cmp1 ≈ 9.8459456377190043e-9 rtol=1e-6
            @test cmp2 ≈ 8.9007105172073563e-11 rtol=1e-6
        end

        @testset "the first zone at the temperature and electron fraction of the reference run" begin
            balance = Radix.ionization_balance(mixture, 100.0, 1e4; radiation=mapped, iterate=false)
            hc = Radix.heating_cooling(mixture, balance, 100.0, 1e4, processes; radiation=mapped)
            # XSTAR's table of the run with lprint=3 (xstar_pow_xi2_eq/README.md): the totals and the terms that are not elements
            @test term(hc, Radix.Compton).heating ≈ 2.37589869e-16 rtol=2e-3
            @test term(hc, Radix.Compton).cooling ≈ 1.95382798e-16 rtol=2e-3
            @test term(hc, Radix.FreeFree).heating ≈ 2.44377767e-26 rtol=1e-3
            @test term(hc, Radix.Bremsstrahlung).cooling ≈ 1.99246730e-16 rtol=1e-3
            @test hc.heating ≈ 5.20663621e-15 rtol=1e-2           # (the database of the tree: 0.1 %, that of the package of XSTAR: 0.1 %)
            @test hc.cooling ≈ 5.77911452e-15 rtol=1e-2
            @test abs(hc.imbalance - (-0.104221974)) < 3e-3
            @test all(e -> e.heating >= 0 && e.cooling >= 0, hc.elements)
            @test hc.heating2 > 0 && hc.cooling2 > 0
            # without the mapping (the flat tail of bremsmap) the Compton terms are a third lower
            plain = Radix.heating_cooling(mixture, balance, 100.0, 1e4, processes; radiation=incident)
            @test term(plain, Radix.Compton).heating ≈ term(hc, Radix.Compton).heating/1.506 rtol=1e-3
        end

        @testset "the iron at log ξ = 1 (a record from the ground level to a level of the next ion)" begin
            # test/reference/xstar_fe_balance: T = 4.5513, x_e = 1. The records of rate type 1 of Fe VII and Fe VIII (`ParPhotoIonize3`, to the levels 3 and 4 of the next ion) used to be left out
            # of the matrix, which put 4-7% of the iron in the wrong ions and made the temperature of the equilibrium 3% low (8% with the database of the tree)
            dir = joinpath(@__DIR__, "reference", "xstar_fe_balance")
            cont = fits(joinpath(dir, "xout_cont1.fits"))[3].data
            ab = fits(joinpath(dir, "xout_abund1.fits"))[2].data
            rad = Radix.map_spectrum(Radix.point_source(Float64.(cont.energy), Float64.(cont.incident)*1e38, 3.16228e20))
            T = 4.5513
            balance = Radix.ionization_balance(mixture, T, 1e4; radiation=rad, xee=1.0, iterate=false)
            iron = findfirst(==(26), mixture.Z)
            for (ion, name) in ((7, :fe_vii), (8, :fe_viii), (9, :fe_ix), (10, :fe_x), (11, :fe_xi), (12, :fe_xii))
                @test balance.fractions[iron][ion] ≈ getproperty(ab, name)[1] rtol=1e-3
            end
            hc = Radix.heating_cooling(mixture, balance, T, 1e4, processes; radiation=rad)
            @test hc.elements[iron].heating ≈ 1.07571729e-15 rtol=1e-3
            @test hc.elements[iron].cooling ≈ 9.67042360e-16 rtol=1e-3
            @test hc.heating ≈ 1.18243531e-14 rtol=1e-3
            @test term(hc, Radix.Compton).heating ≈ 2.37589869e-17 rtol=1e-4
            @test term(hc, Radix.Bremsstrahlung).cooling ≈ 4.14848716e-17 rtol=1e-4
        end

        @testset "thermal equilibrium of the first zone" begin
            eq = fits(joinpath(@__DIR__, "reference", "xstar_pow_xi2_eq", "xout_abund1.fits"))[2].data
            res = Radix.thermal_equilibrium(mixture, 1e4, processes; radiation=mapped, T=100.0)
            @test res.converged && abs(res.imbalance) <= 1e-4
            @test res.xee ≈ eq.x_e[2] rtol=1e-3
            @test res.T ≈ eq.temperature[2] rtol=6e-2             # (2e-4 with the database of the package of XSTAR)
            @test res.nₑ == 1e4*res.xee
            @test res.evaluations < 25
        end
    end
end

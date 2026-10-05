# The transfer through a slab of zones: escape probabilities, line-centre and edge opacities, optical depths, marching.

using FITSFiles

function toy_transfer_tests()
    @testset "Transfer" begin
        @testset "escape probabilities" begin
            # the real pescl and pescv of XSTAR (test/reference/ucalc/drvpesc.f90)
            reference = [(0.0, 0.5, 0.5), (1e-7, 0.5, 4.9999995000000252e-01), (1e-5, 4.9999500003339298e-01, 4.9999500002499991e-01),
                         (1.001e-5, 4.9999499503221262e-01, 4.9999499502504996e-01), (1e-3, 4.9950033316673026e-01, 4.9950024991668751e-01),
                         (0.1, 4.5317311730504545e-01, 4.5241870901797976e-01), (0.5, 3.1606027941427883e-01, 3.0326532985631671e-01),
                         (0.999999, 2.1616632768946523e-01, 1.8393990452553372e-01), (1.0, 2.3507898053286461e-01, 1.8393972058572117e-01),
                         (2.718281828, 6.1045705702968939e-02, 3.2994017937802017e-02), (10.0, 1.4402601402409445e-02, 2.2699964881242427e-05),
                         (1e3, 1.1278741806394859e-04, 4.9999999800209860e-13), (1e5, 1.3772360005478127e-06, 4.9999999800209860e-13),
                         (1e8, 2.3465975655243779e-09, 4.9999999800209860e-13)]
            for (tau, l, v) in reference
                @test Radix.pescl(tau) ≈ l rtol=1e-13
                @test Radix.pescv(tau) ≈ v rtol=1e-13
            end
            @test Radix.pescl(0.0) == 0.5 && Radix.pescv(0.0) == 0.5
            @test issorted([Radix.pescv(t) for t in 0:0.5:30]; rev=true)
            @test Radix.pescl(1.0) > Radix.pescl(1.0 + 1e-9) > Radix.pescl(10.0)        # decreasing across the change of the formula
            @test Radix.pescv(100.0) == 0.5*Float64(1f-12)                              # the floor
        end

        f32 = Float32
        lv(level, E, g) = Radix.AtomicLevel(Int32(13), "", Int32[1, 2, 0, 1, level, 5], f32[E, g, 1, 13.6])
        levels = Radix.levels([lv(1, 0.0, 2), lv(2, 10.2, 8), lv(3, 13.6, 1)]; masses=Dict(5 => 1.0))
        line = Radix.AtomicLine2(Int32(4), "", Int32[2, 1, 1, 5], f32[1215.67, 0.4162, 6.265e8], levels)
        other = Radix.AtomicLine2(Int32(3), "", Int32[2, 1, 1, 5], f32[1215.67, 0.4162, 6.265e8], levels)   # rate type 3: no depth
        elements = Radix.Elements([[line, other]], levels, [5])
        mixture = Radix.Mixture(levels, [1], [0.5], [elements])
        ntot = 1e6
        balance = (; nₑ=1e4, nₕ=0.0, populations=[[0.7, 0.3, 0.0]])

        @testset "opacities and optical depths" begin
            cell = Radix.Cell(1.0, 0.0, 1e4, ntot)
            opacities = Radix.record_opacities(mixture, balance, 1.0, ntot)
            # a line: the cross section at its centre times the density of its lower level (population × abundance × ntot)
            @test opacities[1][1] ≈ Radix.rate(line, cell).opacity*0.7*0.5*ntot
            @test opacities[1][1] > 0 && opacities[1][2] == 0                    # (only the lines of rate type 4 and 9)
            depths = Radix.OpticalDepths(mixture)
            @test depths.inward == [[0.0, 0.0]] && depths.outward == [[0.0, 0.0]]
            Radix.add_zone!(depths, opacities, 1e13)
            @test depths.inward[1][1] ≈ opacities[1][1]*1e13 && depths.outward[1][1] == 0
            Radix.add_zone!(depths, opacities, 2e13; direction=:outward)
            @test depths.outward[1][1] ≈ opacities[1][1]*2e13
            Radix.add_zone!(depths, opacities, 1e13)
            @test depths.inward[1][1] ≈ 2*opacities[1][1]*1e13                    # the zones add
            @test_throws ArgumentError Radix.add_zone!(depths, opacities, 1.0; direction=:sideways)
            # the escape probabilities from the depths
            τ1, τ2 = depths.inward[1][1], depths.outward[1][1]
            escape = Radix.escape_probabilities(mixture, depths)
            @test escape[1][1] == (Radix.pescl(τ1), Radix.pescl(τ2)) && escape[1][2] == (0.5, 0.5)
            cf = 0.3
            @test Radix.escape_probabilities(mixture, depths; cfrac=cf)[1][1] ==
                (Radix.pescl(τ1)*(1 - cf), Radix.pescl(τ2)*(1 - cf) + 2*Radix.pescl(τ1 + τ2)*cf)
            @test Radix.escape_probabilities(mixture, Radix.OpticalDepths(mixture))[1][1] == (0.5, 0.5)    # a thin slab
        end

        @testset "the continuum opacity" begin
            # without photoionization records: Thomson scattering (everywhere) and free-free (in the total only)
            E = Radix.xstar_energy_grid(999)
            rad = Radix.Radiation(E, 1e10 ./ E)
            K = Radix.constants()
            T = 1.0
            continuum = Radix.Continuum(mixture, (Radix.Thomson(), Radix.FreeFree()))
            op = Radix.opacity(continuum, balance, T, ntot, rad, [[0.0, 0.0]])
            thomson = balance.nₑ*K.sigma_thomson
            @test all(op.continuum .== thomson)
            @test op.total ≈ thomson .+ Radix.opacity(Radix.FreeFree(), E, T, balance.nₑ)
            @test op.edges == 0
            @test Radix.opacity(Radix.Continuum(mixture, (Radix.Thomson(0.25), Radix.FreeFree())), balance, T, ntot, rad, [[0.0, 0.0]]).continuum ≈ fill(0.75*thomson, length(E))
            # the free-free opacity is that of the heating
            ff = Radix.opacity(Radix.FreeFree(), E, T, 1e4)
            @test Radix.heating(Radix.FreeFree(), rad, T, 1e4) ≈ K.ergsev*sum((rad.F[k]*ff[k] + rad.F[k - 1]*ff[k - 1])*(E[k] - E[k - 1])/2 for k in 2:length(E))
            @test all(>=(0), ff) && ff[end] > ff[1]*1e-30
        end

        @testset "the lines in the continuum opacity" begin
            cell = Radix.Cell(1.0, 0.0, 1e4, ntot)
            T = 1.0
            edges = Radix.record_opacities(mixture, balance, T, ntot; vturb=10.0)
            # the emissivity of a line: the energy radiated in the decay of the upper level, split by the escape probabilities
            emissivities = Radix.line_emissivities(mixture, balance, T, ntot)
            r = Radix.rate(line, cell; pesc=1.0)
            @test sum(emissivities[1][1]) ≈ 0.3*0.5*ntot*r.fenergy && emissivities[1][1][1] == emissivities[1][1][2]
            @test emissivities[1][2] == (0.0, 0.0)                       # (only the records of type 4)
            thick = Radix.line_emissivities(mixture, balance, T, ntot; escape=[[(0.1, 0.3), (0.5, 0.5)]])
            @test thick[1][1][1]/thick[1][1][2] ≈ 1/3 && sum(thick[1][1]) ≈ 0.3*0.5*ntot*Radix.rate(line, cell; pesc=0.4).fenergy
            @test Radix.line_emissivities(mixture, (; balance..., populations=[[1.0, 0.0, 0.0]]), T, ntot)[1][1] == (0.0, 0.0)
            # the line is added to the bins of the total opacity, and not to the continuum
            E = Radix.xstar_energy_grid(999)
            rad = Radix.Radiation(E, 1e10 ./ E)
            continuum = Radix.Continuum(mixture, (Radix.Thomson(), Radix.FreeFree()); vturb=10.0)
            without = Radix.opacity(continuum, balance, T, ntot, rad, edges)
            with = Radix.opacity(continuum, balance, T, ntot, rad, edges; emissivities)
            added = with.total .- without.total
            @test with.continuum == without.continuum && all(>=(0), added)
            K = Radix.constants()
            e0 = K.hc_eVÅ_single/1215.67
            @test argmax(added) == Radix.nbin(rad, e0) || argmax(added) == Radix.nbin(rad, e0) + 1
            centre = Radix.rate(line, cell; vturb=10.0).opacity*0.7*0.5*ntot
            width = sqrt(((K.thermal_speed/Radix.cm_per_km)*sqrt(T/1.0)*e0/(K.light_speed/Radix.cm_per_km))^2 + (e0*10.0/(K.light_speed/Radix.cm_per_km))^2)
            @test sum(added[k]*(E[k + 1] - E[k]) for k in 1:length(E) - 1) ≈ centre*width rtol=0.05
            # a line weaker than the limit is not added
            weak = (; balance..., populations=[[1e-40, 0.0, 0.0]])
            @test Radix.opacity(continuum, weak, T, ntot, rad, Radix.record_opacities(mixture, weak, T, ntot); emissivities=Radix.line_emissivities(mixture, weak, T, ntot)).total == without.total
        end

        @testset "the transmitted spectrum" begin
            # `heatt`: the spectrum of the source is attenuated and the emission added zone by zone
            K = Radix.constants()
            E = [1.0, 2.0, 3.0, 4.0]
            L = [10.0, 20.0, 30.0, 40.0]
            r, Δr, cfrac = 2e13, 1e12, 0.25
            rad = Radix.point_source(E, L, r)
            total = [1e-14, 1e-12, 5e-12, 1e-30]
            continuum = (; total, continuum=total ./ 2, emissivity=[1e-20 2e-20 3e-20 4e-20; 4e-20 3e-20 2e-20 1e-20], bremsstrahlung=[1e-18, 2e-18, 3e-18, 4e-18])
            start = Radix.Transmitted(L)
            @test start.L == L && all(iszero, start.inward) && all(iszero, start.outward) && all(iszero, start.inward_continuum) && all(iszero, start.outward_continuum)
            after = Radix.transmit(start, rad, continuum, r, Δr; cfrac)
            fpr2 = K.fourpi*r^2
            for k in 1:4
                τ = max(1e-49, total[k])*Δr
                fac = τ > 0.01 ? (1 - exp(-τ))/τ : 1.0
                ε₁ = continuum.emissivity[1, k] + continuum.bremsstrahlung[k]*(1 - cfrac)/2/K.fourpi
                ε₂ = continuum.emissivity[2, k] + continuum.bremsstrahlung[k]*(1 + cfrac)/2/K.fourpi
                @test after.L[k] ≈ max(0.0, L[k] - (rad.F[k]*max(1e-49, total[k]) - K.fourpi*(ε₁ + ε₂))*fac*Δr*fpr2)
                @test after.inward[k] ≈ K.fourpi*ε₁*fac*Δr*fpr2 && after.outward[k] ≈ K.fourpi*ε₂*fac*Δr*fpr2
                τc = max(1e-49, total[k]/2)*Δr
                facc = τc > 0.01 ? (1 - exp(-τc))/τc : 1.0
                @test after.inward_continuum[k] ≈ K.fourpi*ε₁*facc*Δr*fpr2 && after.outward_continuum[k] ≈ K.fourpi*ε₂*facc*Δr*fpr2
            end
            # a thick bin without emission: the flux is what exp(-τ) leaves
            absorbing = (; total=fill(1e-12, 4), continuum=fill(1e-12, 4), emissivity=zeros(2, 4), bremsstrahlung=zeros(4))
            @test Radix.transmit(start, rad, absorbing, r, Δr).L ≈ L .- rad.F .* (1 .- exp.(-1e-12*Δr)) .* fpr2 rtol=1e-12
            # the spectrum does not go below 0 and one zone after another add up
            dark = Radix.transmit(Radix.Transmitted([1.0]), Radix.Radiation([1.0], [1.0]), (; total=[1.0], continuum=[1.0], emissivity=zeros(2, 1), bremsstrahlung=[0.0]), 1e13, 1e13)
            @test dark.L == [0.0]
            twice = Radix.transmit(after, rad, continuum, r, Δr; cfrac)
            @test twice.inward ≈ after.inward .+ K.fourpi .* (continuum.emissivity[1, :] .+ continuum.bremsstrahlung .* (1 - cfrac)/2/K.fourpi) .* [t > 0.01 ? (1 - exp(-t))/t : 1.0 for t in max.(1e-49, total) .* Δr] .* Δr .* fpr2
        end

        @testset "two-photon continua" begin
            levelsx = Radix.levels([lv(1, 0.0, 2), lv(2, 10.2, 8), lv(3, 13.6, 1)]; masses=Dict(5 => 1.0))
            decay = Radix.TwoPhotonDecay(Int32(9), "", Int32[2, 1, 1, 5], f32[8.23], levelsx)
            E = Radix.xstar_energy_grid(999)
            rad = Radix.Radiation(E, zeros(length(E)))
            opacity = Radix.Opacity(length(E))
            abund2 = 1e-3
            Radix.add_two_photon!(opacity, rad, decay, abund2, (0.5, 0.5))
            emax = 10.2
            nbmx = Radix.nbin(rad, emax)
            @test opacity.emissivity[1, :] == opacity.emissivity[2, :]
            @test all(==(0), opacity.emissivity[:, nbmx + 1:end]) && opacity.emissivity[1, 1] == 0
            # the energy of the decay, `A hν` per atom, in all directions: the shape E²(hν - E) normalized
            total = sum(4π*(opacity.emissivity[1, k] + opacity.emissivity[2, k])*(E[k] - E[k - 1]) for k in 2:nbmx)
            @test total ≈ abund2*8.23*emax*(1 + 0.0) rtol=0.1
            shape(k) = E[k]^2*max(0.0, E[nbmx] - E[k])
            @test opacity.emissivity[1, 100]/opacity.emissivity[1, 200] ≈ shape(100)/shape(200)
            # the records with a level out of the table or a rate below 0.01 add nothing
            none = Radix.Opacity(length(E))
            Radix.add_two_photon!(none, rad, Radix.TwoPhotonDecay(Int32(9), "", Int32[2, 1, 1, 5], f32[0.005], levelsx), abund2, (0.5, 0.5))
            Radix.add_two_photon!(none, rad, Radix.TwoPhotonDecay(Int32(9), "", Int32[3, 1, 1, 5], f32[8.23], levelsx), abund2, (0.5, 0.5))
            @test all(==(0), none.emissivity)
        end

        @testset "the thickness of the zones" begin
            # `step`: the smallest of r/steps, the column left and emult/κ of the bins above ectt, not too deep before, with spectrum
            E = [0.5, 5.0, 50.0, 500.0, 5000.0]
            κ = [1e-10, 2e-17, 1e-17, 4e-18, 1e-30]
            L = fill(1e30, 5)
            dpthc = zeros(5)
            args(; r=1e21, depth=0.0, ntot=1e4, column=1e21, emult=1.0, taumax=5.0, ectt=1.0, steps=3) = (r, depth, ntot, column, emult, taumax, ectt, steps)
            thickness(; κ=κ, L=L, dpthc=dpthc, kw...) = Radix.step_thickness(κ, E, L, dpthc, args(; kw...)...)
            @test thickness() == 1/2e-17                                 # (the bin below ectt, with κ = 10⁻¹⁰, is not used)
            @test thickness(; emult=0.5) == 0.5/2e-17
            @test thickness(; ectt=0.1) == 1/1e-10
            @test thickness(; dpthc=[0.0, 5.5, 0.0, 0.0, 0.0]) == 1/1e-17          # the deeper than taumax is not
            @test thickness(; dpthc=[0.0, 5.0, 0.0, 0.0, 0.0]) == 1/2e-17          # (exactly taumax is)
            @test thickness(; L=[1e30, 1e-10, 1e30, 1e30, 1e30]) == 1/1e-17        # spectrum below 10⁻¹² of 10³⁸
            @test thickness(; r=3e16) == 1e16                                    # r / steps
            @test thickness(; depth=9.99e20) ≈ 1e17*1e-3 rtol=1e-6              # the column left
            @test thickness(; κ=zeros(5)) == 1e17                                  # a transparent slab: all the column
        end

        @testset "a slab" begin
            hlevels = Radix.levels([lv(1, 0.0, 2), lv(2, 13.6, 1)]; masses=Dict(5 => 1.0))
            hrate = RecombiningToyRate(Int8(3), Int32(5), 2e-4)
            hmixture = Radix.Mixture(hlevels, [1], [1.0], [Radix.Elements([[hrate]], hlevels, [5])])
            processes = Radix.standard_processes(Radix.Compton([1e-4, 1.0], [1e-3, 1.0], [1.0 2.0; 3.0 4.0]))
            E = Radix.xstar_energy_grid(999)
            L = 1e30 ./ E
            slab = Radix.march_slab(hmixture, 1e4, processes, E, L; r=1e13, column=1e21, T=1.0, iterate=false, diffuse=false)
            # a thin gas (Thomson scattering at about 10⁻²⁰ cm⁻¹): zones of the largest size, r/steps, the first of no thickness, the last ends at the column
            Δ = [z.Δr for z in slab.zones]
            @test Δ[1] == 0 && Δ[2] == 1e13/2                                      # (steps = 2)
            @test slab.zones[1].r == 1e13 && slab.zones[2].r == 1e13 && slab.zones[3].r == 1e13 + Δ[2]
            @test sum(Δ) ≈ 1e17 && all(>=(0), Δ)
            @test sum(Δ[1:end - 1])*1e4 < 1e21 <= sum(Δ)*1e4*(1 + 1e-12)           # the zones go on until the column is reached
            few = Radix.march_slab(hmixture, 1e4, processes, E, L; r=1e13, column=1e21, steps=10, T=1.0, iterate=false, diffuse=false)
            @test few.zones[2].Δr == 1e12 && length(few.zones) > length(slab.zones)
            thin = Radix.march_slab(hmixture, 1e4, processes, E, L; r=1e30, column=1e17, T=1.0, iterate=false, diffuse=false)
            @test [z.Δr for z in thin.zones] == [0.0, 1e13]                          # a thin slab: the column left
        end

        @testset "the source" begin
            E = Radix.xstar_energy_grid(999)
            L = 1e46
            for α in (-1.0, -2.0, 0.5)
                spectrum = Radix.power_law(E, α, L)
                @test spectrum[500]/spectrum[300] ≈ (E[500]/E[300])^α rtol=1e-12
                # the luminosity between 13.6 eV and 13.6 keV is the one asked for
                inside = findall(e -> Float64(13.6f0) <= e <= 1.36e4, E)
                @test sum((spectrum[i] + spectrum[i - 1])*(E[i] - E[i - 1])/2 for i in inside)*Radix.constants().ergsev ≈ L rtol=1e-12
            end
            @test Radix.power_law(E, -1.0, 2L) ≈ 2 .* Radix.power_law(E, -1.0, L) rtol=1e-12
            @test Radix.power_law([0.005, 0.5, 5.0, 50.0, 500.0], -1.0, L)[1] < 1e-20*Radix.power_law([0.005, 0.5, 5.0, 50.0, 500.0], -1.0, L)[2]   # (the floor below 0.01 eV)
            @test Radix.source_distance(1e46, 10.0, 1e4) ≈ 3.1622776601683794e20 && Radix.source_distance(4e46, 10.0, 1e4) ≈ 2*Radix.source_distance(1e46, 10.0, 1e4)
        end

        @testset "the escape of one element" begin
            per_element = [[(0.1, 0.2), (0.3, 0.4)]]
            @test Radix.element_escape(per_element, 1) == [(0.1, 0.2), (0.3, 0.4)]
            @test Radix.element_escape([(0.1, 0.2), (0.3, 0.4)], 1) == [(0.1, 0.2), (0.3, 0.4)]    # one vector of pairs for all
            @test Radix.element_escape(nothing, 3) === nothing
            f = coef -> (0.5, 0.5)
            @test Radix.element_escape(f, 2) === f
        end

        @testset "marching" begin
            # one ion, ionization at 1 s⁻¹ and recombination at α nₑ: nothing for the transfer to do (no lines): the zones are solved
            # in turn at their own radius, with the temperature and the electron fraction of the one before
            hlevels = Radix.levels([lv(1, 0.0, 2), lv(2, 13.6, 1)]; masses=Dict(5 => 1.0))
            hrate = RecombiningToyRate(Int8(3), Int32(5), 2e-4)
            helements = Radix.Elements([[hrate]], hlevels, [5])
            hmixture = Radix.Mixture(hlevels, [1], [1.0], [helements])
            processes = Radix.standard_processes(Radix.Compton([1e-4, 1.0], [1e-3, 1.0], [1.0 2.0; 3.0 4.0]))
            E = Radix.xstar_energy_grid(999)
            L = 1e30 ./ E
            run = Radix.march_zones(hmixture, 1e4, processes, E, L, [(1e13, 1e12), (2e13, 3e12)]; T=1.0, iterate=false, diffuse=false)
            @test length(run.zones) == 2
            @test [z.r for z in run.zones] == [1e13, 2e13] && [z.Δr for z in run.zones] == [1e12, 3e12]
            @test run.zones[1].radiation.F == Radix.map_spectrum(Radix.point_source(E, L, 1e13)).F
            @test run.zones[2].radiation.F[10] ≈ run.zones[1].radiation.F[10]/4                  # 1/r²
            @test all(z -> z.T == 1.0 && z.xee == 1.0, run.zones)
            @test run.depths.inward == [[0.0]] && run.zones[1].escape == [[(0.5, 0.5)]]
            # the continuum is attenuated by the depth of the zones before (here Thomson scattering): exp(-Σ opacity Δr)
            thick = Radix.march_zones(hmixture, 1e4, processes, E, L, [(1e13, 1e17), (1e13, 1e17)]; T=1.0, iterate=false, diffuse=false)
            thomson = thick.zones[1].nₑ*Radix.constants().sigma_thomson
            @test thick.zones[1].opacity.total[10] ≈ thomson + Radix.opacity(Radix.FreeFree(), E, 1.0, thick.zones[1].nₑ)[10]
            @test thick.zones[2].radiation.F[10] ≈ thick.zones[1].radiation.F[10]*exp(-thick.zones[1].opacity.total[10]*1e17) rtol=1e-6      # (`heatt` steps linearly in a thin zone)
            @test thick.dpthc[10] ≈ 2*thick.zones[1].opacity.total[10]*1e17
            @test thick.dpthcont[10] ≈ 2*thomson*1e17                                      # (free-free is not in the continuum)
            @test thick.dpthc[10] > thick.dpthcont[10]
            @test Radix.march_zones(hmixture, 1e4, processes, E, L, [(1e13, 1e17), (1e13, 1e17)]; T=1.0, iterate=false, diffuse=false, lines=false).dpthc == thick.dpthc    # (no lines here)
            # the diffuse emission (here the bremsstrahlung of the zone) adds to the spectrum that the next zone has
            emitting = Radix.march_zones(hmixture, 1e4, processes, E, L, [(1e13, 1e17), (1e13, 1e17)]; T=1.0, iterate=false)
            @test emitting.zones[2].radiation.F[10] > thick.zones[2].radiation.F[10]
            @test emitting.spectrum.inward_continuum[10] > 0 && emitting.spectrum.L[10] > thick.spectrum.L[10]
            flat = Radix.march_zones(hmixture, 1e4, processes, E, L, [(1e13, 1e17), (1e13, 1e17)]; T=1.0, iterate=false, attenuate=false)
            @test flat.zones[2].radiation.F == flat.zones[1].radiation.F
            @test run.zones[1].fractions[1][2] ≈ 1/3                                      # x₂ = 1/(1 + α nₑ) with α nₑ = 2
        end
    end
end

# the reference run of XSTAR (test/reference/xstar_pow_xi2): two zones of 5×10¹² cm at 10¹³ and 1.5×10¹³ cm, the second the one of the
# table of ion fractions. Its files give the optical depths of the lines (depth_inward) and of the recombination edges.
function transfer_balance_tests(db)
    @testset "Transfer of the database" begin
        levels = Radix.levels(db)
        coheat = joinpath(dirname(get(ENV, "RADIX_ATDB", "")), "coheat.dat")
        isfile(coheat) || (@info "Skipping the transfer tests (coheat.dat next to atdb.fits is needed)"; return)
        processes = Radix.standard_processes(Radix.load(Radix.Compton, coheat))
        mixture = Radix.Mixture(db, levels; multiplier=Dict(3 => 0.0, 4 => 0.0, 5 => 0.0))
        dir = joinpath(@__DIR__, "reference", "xstar_pow_xi2")
        spectrum = fits(joinpath(dir, "xout_cont1.fits"))[3].data
        lines = fits(joinpath(dir, "xout_lines1.fits"))[3].data
        edges = fits(joinpath(dir, "xout_rrc1.fits"))[3].data
        abundances = fits(joinpath(dir, "xout_abund1.fits"))[2].data
        E, L = Float64.(spectrum.energy), Float64.(spectrum.incident)*1e38
        run = Radix.march_zones(mixture, 1e4, processes, E, L, [(1e13, 5e12), (1.5e13, 5e12)]; T=100.0, iterate=false)
        ionlabel = Dict(Int(r.ion) => strip(r.label) for r in db if r isa Radix.Ion)
        label(ion, level) = strip(levels[(Int(ion), Int(level))].label)

        @testset "the two zones against the table of the reference run" begin
            for (z, row) in zip(run.zones, (2, 3))
                @test abs(z.imbalance - abundances.frac_heat_error[row]) < 2e-3
            end
            # the optical depths are those of the first zone's lines, which are tiny: the second zone is as without trapping
            @test maximum(maximum, run.depths.inward) < 6e-3
            @test run.zones[2].heating ≈ Radix.heating_cooling(mixture, Radix.ionization_balance(mixture, 100.0, 1e4;
                radiation=run.zones[2].radiation, iterate=false), 100.0, 1e4, processes; radiation=run.zones[2].radiation).heating rtol=1e-5
        end

        @testset "optical depths of the lines" begin
            depth = Dict{Tuple{String, String, String}, Float64}()
            for (k, el) in enumerate(mixture.elements), j in eachindex(el.rates)
                c = el.rates[j]
                (c isa Radix.AtomicLine2 && c.rtype in Radix.line_rate_types) || continue
                lo, up = el.lo[j], el.up[j]
                lo == 0 && continue
                key = (ionlabel[Int(c.ion)], label(c.ion, c.transition.lower), label(c.ion, c.transition.upper))
                run.depths.inward[k][j] > 0 && (depth[key] = get(depth, key, 0.0) + run.depths.inward[k][j])
            end
            ratios = Float64[]
            for i in eachindex(lines.ion)
                ref = Float64(lines.depth_inward[i])
                ref > 1e-9 || continue
                key = (strip(lines.ion[i]), strip(lines.lower_level[i]), strip(lines.upper_level[i]))
                mine = get(depth, key, nothing)
                mine === nothing && (mine = get(depth, (key[1], key[3], key[2]), 0.0))      # (the levels in the other order)
                push!(ratios, mine/ref)
            end
            @test length(ratios) > 400
            @test all(r -> abs(r - 1) < 2e-2, ratios)
            @test abs(sort(ratios)[length(ratios) ÷ 2] - 1) < 2e-3                           # the median
        end

        @testset "the continuum depth and the emission of a thick slab" begin
            # test/reference/xstar_thick_slab: 10²¹ cm⁻² at log ξ = 1, T = 10⁵ K, without the vturb of XSTAR (see its README). The zones are those of XSTAR: the rows
            # of its table are the states at their inner edges (`delta_r`), the last being the end of the slab
            dir = joinpath(@__DIR__, "reference", "xstar_thick_slab")
            slab = fits(joinpath(dir, "xout_cont1.fits"))[3].data
            Es, Is, Ts = Float64.(slab.energy), Float64.(slab.incident), Float64.(slab.transmitted)
            slabab = fits(joinpath(dir, "xout_abund1.fits"))[2].data
            r0 = 3.16228e20
            depth = Float64.(slabab.delta_r)[2:end]
            zones = [(r0 + depth[k], depth[k + 1] - depth[k]) for k in 1:length(depth) - 1]
            thickrun = Radix.march_zones(mixture, 1e4, processes, Es, Is*1e38, zones; T=10.0, iterate=false, vturb=0.0, luminous=false)
            τx = -log.(Ts ./ Is)
            ratio = thickrun.dpthcont ./ τx
            @test all(τx .> 6e-4)
            @test abs(sort(ratio)[length(ratio) ÷ 2] - 1) < 5e-3                 # the median of the 999 bins
            @test all(r -> abs(r - 1) < 6e-2, ratio)
            i = argmin(abs.(Es .- 13.6)); @test thickrun.dpthcont[i] ≈ τx[i] rtol=2e-3           # Thomson scattering
            i = argmin(abs.(Es .- 54.7)); @test thickrun.dpthcont[i] ≈ τx[i] rtol=4e-2           # the edge of He II
            # the states at the inner edges of the zones follow the table: the first is that of the table, the others differ by the diffuse emission that the zones add to the
            # radiation of those after them (`transmit`), which XSTAR's table has: 0.1% in He II without the emission 2.7%
            @test thickrun.zones[1].fractions[2][2] ≈ slabab.he_ii[2] rtol=1e-4
            @test thickrun.zones[1].fractions[findfirst(==(8), mixture.Z)][8] ≈ slabab.o_viii[2] rtol=1e-4
            # the thickness of the zones of XSTAR's `step` (its log: log N = 20.26, 20.56, 20.74, 20.86, 20.95, 21.00 at the inner edges of the zones after the first two and at the end)
            auto = Radix.slab_model(mixture, processes; density=1e4, column=1e21, logξ=1.0, luminosity=1e8*1e38, α=-1.0, emult=1.0, taumax=5.0, steps=3, T=10.0, iterate=false, vturb=0.0, luminous=false)
            @test length(auto.zones) == 7 && auto.zones[1].Δr == 0
            @test auto.r ≈ r0 rtol=1e-6                              # the source of the reference run: E⁻¹ with the luminosity 10⁴⁶ erg/s
            @test auto.E ≈ Es rtol=1e-5
            @test auto.incident ./ 1e38 ≈ Is rtol=1e-5
            @test abs(sort(auto.transmitted ./ 1e38 ./ Ts)[end ÷ 2] - 1) < 5e-3 && all(abs.(auto.transmitted ./ 1e38 ./ Ts .- 1) .< 6e-2)
            @test [z.r - auto.r for z in auto.zones][1:6] ≈ [0; depth[1:5]] atol=3e-4*1.8e16 rtol=1e-3
            @test round.(log10.(1e4 .* cumsum([z.Δr for z in auto.zones])[2:end]); digits=2) == [20.26, 20.56, 20.74, 20.86, 20.95, 21.0]
            @test sum(z.Δr for z in auto.zones) ≈ 1e17
            for k in 2:length(zones)
                @test thickrun.zones[k].fractions[2][2] ≈ slabab.he_ii[k + 1] rtol=1.5e-3
                @test thickrun.zones[k].fractions[findfirst(==(8), mixture.Z)][8] ≈ slabab.o_viii[k + 1] rtol=3e-4
            end
            # the diffuse emission of the zones, summed, is the emission columns of the transmitted spectrum of XSTAR (the bremsstrahlung, the recombination continua and
            # the two-photon continua): the bins with at least a thousandth of the strongest
            for (mine, reference) in ((thickrun.spectrum.inward_continuum, slab.emit_inward), (thickrun.spectrum.outward_continuum, slab.emit_outward))
                reference = Float64.(reference)
                strong = findall(>(1e-3*maximum(reference)), reference)
                bins = (mine[strong] ./ 1e38) ./ reference[strong]
                @test length(strong) > 500
                @test abs(sort(bins)[length(bins) ÷ 2] - 1) < 3e-3
                @test all(r -> abs(r - 1) < 6e-2, bins)
                @test sum(mine)/1e38 ≈ sum(reference) rtol=1e-2
            end
            # the end of the spectrum is the transmitted spectrum, where the lines do not absorb the bins that the edges do
            # the XSTAR of the reference does not absorb the continuum at all when vturb > 0 (gsmooth2)
            smoothed = fits(joinpath(dir, "xout_cont1_vturb1.fits"))[3].data
            @test count(==(1), Float64.(smoothed.transmitted) ./ Float64.(smoothed.incident)) > 700
        end

        @testset "the luminosities of the lines and the recombination edges of the thick slab" begin
            # the zones of `slab_model` (the seven of XSTAR) with its table of the lines (`xout_lines1.fits`: emit_inward + emit_outward in 10³⁸ erg/s, for the lines above 10²) and the total of
            # its table of the recombination edges (`xout_rrc1.fits`, 6 MB, not kept): 807150.6. The run of the table has the hydrogen of the database of the package of XSTAR
            slab = Radix.slab_model(mixture, processes; density=1e4, column=1e21, logξ=1.0, luminosity=1e8*1e38, emult=1.0, taumax=5.0, steps=3, T=10.0, iterate=false, vturb=0.0)
            ionlabel = Dict(Int(r.ion) => strip(r.label) for r in db if r isa Radix.Ion)
            label(ion, level) = strip(levels[(Int(ion), Int(level))].label)
            emission = Dict{Tuple{String, String, String}, Float64}()
            edge_total = 0.0
            for (k, el) in enumerate(mixture.elements), j in eachindex(el.rates)
                c = el.rates[j]
                total = slab.luminosities.inward[k][j] + slab.luminosities.outward[k][j]
                if c isa Radix.AtomicLine2 && c.rtype == Radix.line_data_type && el.lo[j] != 0 && total > 0
                    key = (ionlabel[Int(c.ion)], label(c.ion, c.transition.lower), label(c.ion, c.transition.upper))
                    emission[key] = get(emission, key, 0.0) + total
                elseif c.rtype == Radix.edge_rate_type && c isa Radix.opacity_edge_types
                    edge_total += total
                end
            end
            @test edge_total ≈ 807150.611744286 rtol=1e-3
            table = fits(joinpath(@__DIR__, "reference", "xstar_thick_slab", "xout_lines1.fits"))[3].data
            ratios = Float64[]
            for i in eachindex(table.ion)
                reference = Float64(table.emit_inward[i]) + Float64(table.emit_outward[i])
                reference > 1e2 || continue
                key = (strip(table.ion[i]), strip(table.lower_level[i]), strip(table.upper_level[i]))
                mine = get(emission, key, get(emission, (key[1], key[3], key[2]), nothing))
                @test mine !== nothing
                mine === nothing && continue
                key[1] == "h_i" && continue                           # (the hydrogen of the database of the tree differs: 2×)
                push!(ratios, mine/reference)
            end
            @test length(ratios) > 500
            @test abs(sort(ratios)[length(ratios) ÷ 2] - 1) < 1e-3
            @test all(r -> abs(r - 1) < 6e-2, ratios)
        end

        @testset "optical depths of the recombination edges" begin
            # the edges of H-like and He-like ions, which have a single record each (the table has no parent level to tell the others)
            wanted = Dict(("o_viii", "1s1.2S_1/2") => 871.4, ("c_vi", "1s1.2S_1/2") => 490.0, ("he_ii", "1s1.2S_1/2") => 54.42,
                          ("si_xiii", "1s2.1S_0") => 2438.0, ("ne_x", "1s1.2S_1/2") => 1362.0, ("mg_xii", "1s1.2S_1/2") => 1963.0)
            for ((ion, level), eth) in wanted
                refs = [Float64(edges.depth_outward[i]) for i in eachindex(edges.ion)
                        if strip(edges.ion[i]) == ion && strip(edges.level[i]) == level && abs(Float64(edges.energy[i]) - eth) < 1]
                @test length(refs) == 1
                mine = 0.0
                for (k, el) in enumerate(mixture.elements), j in eachindex(el.rates)
                    c = el.rates[j]
                    c isa Radix.ParPhotoIonize2 || continue
                    el.lo[j] == 0 && continue
                    (ionlabel[Int(c.ion)] == ion && label(c.ion, c.level) == level) || continue
                    mine += run.depths.inward[k][j]
                end
                @test mine ≈ only(refs) rtol=1e-2
            end
        end
    end
end

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
            slab = Radix.iterate(hmixture, 1e4, processes, E, L; r=1e13, column=1e21, T=1.0, iterate=false, diffuse=false)
            # a thin gas (Thomson scattering at about 10⁻²⁰ cm⁻¹): zones of the largest size, r/steps, the first of no thickness, the last ends at the column
            Δ = [z.Δr for z in slab.zones]
            @test Δ[1] == 0 && Δ[2] == 1e13/2                                      # (steps = 2)
            @test slab.zones[1].r == 1e13 && slab.zones[2].r == 1e13 && slab.zones[3].r == 1e13 + Δ[2]
            @test sum(Δ) ≈ 1e17 && all(>=(0), Δ)
            @test sum(Δ[1:end - 1])*1e4 < 1e21 <= sum(Δ)*1e4*(1 + 1e-12)           # the zones go on until the column is reached
            # another source
            hot = Radix.slab(hmixture, processes; ncn2=999, density=1e4, column=1e17, logξ=2.0, luminosity=1e35, spectrum=(E, L) -> Radix.blackbody(E, 0.5, L), T=1.0, iterate=false, diffuse=false)
            @test hot.incident == Radix.blackbody(hot.E, 0.5, 1e35) && hot.r ≈ sqrt(1e35/(1e2*1e4)) && hot.zones[1].radiation.F[10] ≈ Radix.map_spectrum(Radix.point_source(hot.E, hot.incident, hot.r)).F[10]
            # constant pressure: the density of each zone is P/(k T) at the temperature of the zone (here fixed), and the column adds up with it
            pressure = 1.38e-12*1e4*1.0                                          # n = 10⁴ at T = 1 (10⁴ K)
            @test Radix.ConstantPressure(pressure)(1.0) ≈ 1e4 rtol=1e-6
            @test Radix.ConstantPressure(pressure)(4.0) ≈ 2.5e3 rtol=1e-6
            @test Radix.gas_density(5.0, 2.0) == 5.0
            @test Radix.gas_density(Radix.ConstantPressure(pressure), 2.0) ≈ 5e3 rtol=1e-6
            cp = Radix.iterate(hmixture, Radix.ConstantPressure(pressure), processes, E, L; r=1e13, column=1e21, T=1.0, iterate=false, diffuse=false)
            @test all(z -> isapprox(z.ntot, 1e4; rtol=1e-6) && z.T == 1.0, cp.zones)
            @test isapprox([z.Δr for z in cp.zones], [z.Δr for z in slab.zones]; rtol=1e-5)          # (the same as the constant density of 10⁴ at the same temperature)
            hot = Radix.iterate(hmixture, Radix.ConstantPressure(pressure), processes, E, L; r=1e13, column=1e21, T=2.0, iterate=false, diffuse=false)
            @test all(z -> isapprox(z.ntot, 5e3; rtol=1e-6), hot.zones)
            @test sum(z.Δr for z in hot.zones) ≈ 1e21/5e3 rtol=1e-6        # half the density: twice the length for the column
            @test_throws ArgumentError Radix.slab(hmixture, processes; ncn2=999, column=1e21, logξ=2.0, luminosity=1e35, T=1.0, iterate=false)
            @test Radix.source_distance_pressure(1e46, 10.0, 1e-7) ≈ sqrt(1e46/(Radix.constants().fourpi*Float64(2.99792458f10)*1e-7*10.0))
            # a density that varies as a power of the distance from the source (`radexp`): the density of a zone is that of its inner edge, and the column adds up with it
            d = Radix.PowerLawDensity(1e4, 1e13, -2.0)
            @test Radix.gas_density(d, 5.0) == 1e4
            @test Radix.gas_density(d, 5.0, 2e13) ≈ 2.5e3
            @test Radix.gas_density(1e4, 5.0, 2e13) == 1e4
            radial = Radix.iterate(hmixture, d, processes, E, L; r=1e13, column=5e16, T=1.0, iterate=false, diffuse=false)        # (the column of n ∝ r⁻² from 10¹³ cm adds up to at most 6.7×10¹⁶ with the zones of r/2)
            @test all(z -> isapprox(z.ntot, 1e4*(z.r/1e13)^-2; rtol=1e-12), radial.zones)
            outer(z) = Radix.gas_density(d, 1.0, z.r + z.Δr)
            @test sum(outer(z)*z.Δr for z in radial.zones) ≈ 5e16 rtol=1e-9           # (the column adds up with the density at the outer edge of each zone, as XSTAR does)
            @test length(radial.zones) > 1
            @test_throws ArgumentError Radix.iterate(hmixture, d, processes, E, L; r=1e13, column=6.7e16, T=1.0, iterate=false, diffuse=false)      # (out of reach: the zones are at most r/2 thick and the density falls as r⁻²)
            @test radial.zones[end].ntot < 1e4
            model = Radix.slab(hmixture, processes; ncn2=999, density=1e4, radexp=-2.0, column=1e16, logξ=2.0, luminosity=1e35, T=1.0, iterate=false, diffuse=false)
            @test model.zones[1].ntot ≈ 1e4
            @test model.zones[end].ntot ≈ 1e4*(model.zones[end].r/model.r)^-2 rtol=1e-12
            @test Radix.slab(hmixture, processes; ncn2=999, density=1e4, radexp=0.0, column=1e21, logξ=2.0, luminosity=1e35, T=1.0, iterate=false, diffuse=false).zones[end].ntot == 1e4
            few = Radix.iterate(hmixture, 1e4, processes, E, L; r=1e13, column=1e21, steps=10, T=1.0, iterate=false, diffuse=false)
            @test few.zones[2].Δr == 1e12 && length(few.zones) > length(slab.zones)
            thin = Radix.iterate(hmixture, 1e4, processes, E, L; r=1e30, column=1e17, T=1.0, iterate=false, diffuse=false)
            @test [z.Δr for z in thin.zones] == [0.0, 1e13]                          # a thin slab: the column left
        end

        @testset "the source" begin
            E = Radix.xstar_energy_grid(999)
            L = 1e46
            for α in (-1.0, -2.0, 0.5)
                spectrum = Radix.powerlaw(E, α, L)
                @test spectrum[500]/spectrum[300] ≈ (E[500]/E[300])^α rtol=1e-12
                # the luminosity between 13.6 eV and 13.6 keV is the one asked for
                inside = findall(e -> Float64(13.6f0) <= e <= 1.36e4, E)
                @test sum((spectrum[i] + spectrum[i - 1])*(E[i] - E[i - 1])/2 for i in inside)*Radix.constants().ergsev ≈ L rtol=1e-12
            end
            @test Radix.powerlaw(E, -1.0, 2L) ≈ 2 .* Radix.powerlaw(E, -1.0, L) rtol=1e-12
            @test Radix.powerlaw([0.005, 0.5, 5.0, 50.0, 500.0], -1.0, L)[1] < 1e-20*Radix.powerlaw([0.005, 0.5, 5.0, 50.0, 500.0], -1.0, L)[2]   # (the floor below 0.01 eV)
            # the thermal spectra and the real ones of XSTAR (test/reference/xstar_sources): the bins with at least 10⁻³⁰ of the maximum, to the precision of the files
            sources = joinpath(@__DIR__, "reference", "xstar_sources")
            for (name, f, T) in (("bbody_0.5", Radix.blackbody, 0.5), ("bbody_0.0003", Radix.blackbody, 0.0003), ("brems_0.1", Radix.thermal, 0.1))
                reference = Float64.(fits(joinpath(sources, name * ".fits"))[3].data.incident)
                mine = Radix.constants(Radix.ucalc()) do        # (XSTAR's k: 0.861707 eV per 10⁴ K)
                    f(Radix.xstar_energy_grid(999), T, 1e35) ./ 1e38
                end
                strong = findall(>(1e-30*maximum(reference)), reference)
                @test length(strong) > 800
                @test all(abs.(mine[strong] ./ reference[strong] .- 1) .< 1e-4)
            end
            # a table: a power law is exact in the logarithms, the photons and the logarithms of the units are as they are read, and nothing is outside the table
            table_E = [0.5, 3.0, 40.0, 700.0, 9000.0, 3e5]
            reference = Radix.powerlaw(E, -1.5, L)
            inside = findall(e -> 0.5 <= e <= 3e5, E)
            @test Radix.tabulated(E, table_E, table_E .^ -1.5, L)[inside] ≈ reference[inside] rtol=1e-10
            @test Radix.tabulated(E, table_E, table_E .^ -2.5, L; units=:photons)[inside] ≈ reference[inside] rtol=1e-10
            @test Radix.tabulated(E, table_E, -1.5 .* log10.(table_E), L; units=:log10)[inside] ≈ reference[inside] rtol=1e-10
            @test all(==(0), Radix.tabulated(E, table_E, table_E .^ -1.5, L)[setdiff(eachindex(E), inside)])
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
            run = Radix.iterate(hmixture, 1e4, processes, E, L, [(1e13, 1e12), (2e13, 3e12)]; T=1.0, iterate=false, diffuse=false)
            @test length(run.zones) == 2
            @test [z.r for z in run.zones] == [1e13, 2e13] && [z.Δr for z in run.zones] == [1e12, 3e12]
            @test run.zones[1].radiation.F == Radix.map_spectrum(Radix.point_source(E, L, 1e13)).F
            @test run.zones[2].radiation.F[10] ≈ run.zones[1].radiation.F[10]/4                  # 1/r²
            @test all(z -> z.T == 1.0 && z.xee == 1.0, run.zones)
            @test run.depths.inward == [[0.0]] && run.zones[1].escape == [[(0.5, 0.5)]]
            # the continuum is attenuated by the depth of the zones before (here Thomson scattering): exp(-Σ opacity Δr)
            thick = Radix.iterate(hmixture, 1e4, processes, E, L, [(1e13, 1e17), (1e13, 1e17)]; T=1.0, iterate=false, diffuse=false)
            thomson = thick.zones[1].nₑ*Radix.constants().sigma_thomson
            @test thick.zones[1].opacity.total[10] ≈ thomson + Radix.opacity(Radix.FreeFree(), E, 1.0, thick.zones[1].nₑ)[10]
            @test thick.zones[2].radiation.F[10] ≈ thick.zones[1].radiation.F[10]*exp(-thick.zones[1].opacity.total[10]*1e17) rtol=1e-6      # (`heatt` steps linearly in a thin zone)
            @test thick.dpthc[10] ≈ 2*thick.zones[1].opacity.total[10]*1e17
            @test thick.dpthcont[10] ≈ 2*thomson*1e17                                      # (free-free is not in the continuum)
            @test thick.dpthc[10] > thick.dpthcont[10]
            @test Radix.iterate(hmixture, 1e4, processes, E, L, [(1e13, 1e17), (1e13, 1e17)]; T=1.0, iterate=false, diffuse=false, lines=false).dpthc == thick.dpthc    # (no lines here)
            # the diffuse emission (here the bremsstrahlung of the zone) adds to the spectrum that the next zone has
            emitting = Radix.iterate(hmixture, 1e4, processes, E, L, [(1e13, 1e17), (1e13, 1e17)]; T=1.0, iterate=false)
            @test emitting.zones[2].radiation.F[10] > thick.zones[2].radiation.F[10]
            @test emitting.spectrum.inward_continuum[10] > 0 && emitting.spectrum.L[10] > thick.spectrum.L[10]
            flat = Radix.iterate(hmixture, 1e4, processes, E, L, [(1e13, 1e17), (1e13, 1e17)]; T=1.0, iterate=false, attenuate=false)
            @test flat.zones[2].radiation.F == flat.zones[1].radiation.F
            @test run.zones[1].fractions[1][2] ≈ 1/3                                      # x₂ = 1/(1 + α nₑ) with α nₑ = 2
        @testset "a saturated line in the bins: the profile inside the cell, not the mean opacity" begin
            # a line narrower than a bin (the thermal width is 6×10⁻⁵ of the energy, the bin 1.6×10⁻² at 999 points) absorbs less than its mean opacity says when it is saturated:
            # the transmission of the bin is the average of exp(-τ) over the bin, not exp of the average of τ. `resolved_transmission` gives the first
            K = Radix.constants()
            e0, a, S, Tline = 1500.0, 1e-3, 120.0, 1.0                         # a centre depth of S/√π ≈ 68, a hydrogen-like mass at 10⁴ K
            dele = e0*K.thermal_speed/Radix.cm_per_km*sqrt(Tline)/(K.light_speed/Radix.cm_per_km)
            profile(x) = S*Radix.voigt(x, a)/Radix.line_profile_norm
            totals = Float64[]
            for n in (999, 3000, 10000)
                E = Radix.xstar_energy_grid(n)
                rad = Radix.Radiation(E, zeros(n))
                op = Radix.Opacity(n)
                Radix.add_line!(op, rad, S, 0.0, 0.0, K.hc_eVÅ_single/e0, 0.0, Tline, 1.0, a*dele*K.fourpi; lfast=2)
                absorbers = Radix.LineAbsorbers()
                absorbers.index[(1, 1)] = 1
                push!(absorbers.e0, e0); push!(absorbers.S, S); push!(absorbers.widths, S*dele); push!(absorbers.damping, S*a)
                T = Radix.resolved_transmission(E, absorbers, zeros(n))
                # the average of exp(-τ) over each cell, on a fine grid
                cell(k) = (k == 1 ? E[1] : (E[k - 1] + E[k])/2, k == n ? E[n] : (E[k] + E[k + 1])/2)
                k0 = argmin(abs.(E .- e0))
                exact = map(k0 - 3:k0 + 3) do k
                    lo, hi = cell(k)
                    xs = range(lo, hi; length=40001)
                    sum((exp(-profile((xs[i] - e0)/dele)) + exp(-profile((xs[i + 1] - e0)/dele)))/2 for i in 1:40000)/40000
                end
                @test T[k0 - 3:k0 + 3] ≈ exact atol=5e-3
                weights = [((k == n ? E[n] : (E[k] + E[k + 1])/2) - (k == 1 ? E[1] : (E[k - 1] + E[k])/2)) for k in 1:n]
                push!(totals, sum((1 .- T) .* weights))
                n == 999 && @test abs(sum((1 .- exp.(-op.total)) .* weights)/totals[end] - 1) > 0.1           # (exp of the mean opacity of the bin is not the mean of exp)
            end
            @test maximum(totals)/minimum(totals) < 1.02                         # the equivalent width does not depend on the bins
        end

        @testset "an edge in the bins of the continuum: the average of each cell, not the value at the grid points" begin
            # the recombination continuum of a K edge falls in about kT (9 eV at 10⁵ K) and is one or two bins wide on a grid of 999 points, so XSTAR's value at the grid points
            # (and its first bin at the nearest one to the threshold) gives an integral that depends on where the threshold falls in a bin (a factor 4 at 999 bins for the Ne X edge): `exact_binning`
            K = Radix.constants()
            ε = [0.0, 0.01, 0.05, 0.2, 1.0, 5.0, 30.0]                                  # Ry above the threshold
            σ = 6e-19 ./ (1 .+ 8 .* ε).^3
            eth, T, rnist, ptmp, abund, ntot = 1362.2, 10.0, 2e-3, (0.4, 0.6), (3e-5, 2e-5), 1e4
            kT = K.kB_cgs*T*Radix.T_unit/K.ergsev
            nodes = eth .+ ε .* K.Ry_eV
            # the integral of the table and of the emission integrand on a very fine grid (the table is linear between its points)
            fine(f) = sum(begin
                    xs = range(nodes[j], min(nodes[j + 1], eth + 60kT); length=20001)
                    nodes[j] >= eth + 60kT ? 0.0 : sum((f(xs[i], j) + f(xs[i + 1], j))/2*step(xs) for i in 1:20000)
                end for j in 1:length(ε) - 1)
            table(x, j) = σ[j] + (σ[j + 1] - σ[j])*(x - nodes[j])/(nodes[j + 1] - nodes[j])
            absorption = sum((σ[j] + σ[j + 1])/2*(nodes[j + 1] - nodes[j]) for j in 1:length(ε) - 1)
            emission = abund[2]*ntot*rnist*K.bb_coeff*2*fine((x, j) -> min(Radix.bb_energy_cap, x)^3*table(x, j)*exp(-(x - eth)/kT))
            for n in (999, 2000, 3000, 10000)
                E = Radix.xstar_energy_grid(n)
                op = Radix.Opacity(n)
                Radix.exact_edge_arrays!(op, E, eth, ε, σ, K, kT, rnist, ptmp, abund[1], abund[2], ntot, 2)
                weights = [((k == n ? E[n] : (E[k] + E[k + 1])/2) - (k == 1 ? E[1] : (E[k - 1] + E[k])/2)) for k in 1:n]
                @test sum(op.total .* weights) ≈ abund[1]*ntot*absorption rtol=1e-9
                @test op.continuum == op.total
                @test sum((op.emissivity[1, :] .+ op.emissivity[2, :]) .* weights) ≈ emission rtol=2e-4
                @test sum(op.emissivity[2, :] .* weights)/sum(op.emissivity[1, :] .* weights) ≈ ptmp[2]/ptmp[1] rtol=1e-12
                below = [k for k in 1:n - 1 if (E[k] + E[k + 1])/2 <= eth]                  # (the cells under the threshold have nothing, the one with it has the part above)
                @test all(op.total[below] .== 0) && all(op.emissivity[:, below] .== 0) && all(op.emissivity .>= 0)
            end
            # XSTAR's rule does not conserve it: the integral of its values at the grid points changes with the grid by tens of percent
            old = Radix.constants(Radix.Constants(K; exact_binning=false)) do
                map((999, 2000, 3000)) do n
                    E = Radix.xstar_energy_grid(n)
                    op = Radix.Opacity(n)
                    Radix.photoionization_integrals(Radix.Radiation(E, ones(n)), eth, ε, σ, T, rnist, ptmp; abund, ntot, lfast=2, opacity=op)
                    sum(((op.emissivity[1, k] + op.emissivity[2, k] + op.emissivity[1, k + 1] + op.emissivity[2, k + 1])/2)*(E[k + 1] - E[k]) for k in 1:n - 1)
                end
            end
            @test maximum(old)/minimum(old) > 1.1
        end
        end
    end
end

# the reference run of XSTAR (test/reference/xstar_pow_xi2): two zones of 5×10¹² cm at 10¹³ and 1.5×10¹³ cm, the second the one of the
# table of ion fractions. Its files give the optical depths of the lines (depth_inward) and of the recombination edges.
function transfer_balance_tests(db, known=nothing)
    @testset verbose=true "Transfer of the database" begin
        levels = Radix.levels(db)
        coheat = joinpath(dirname(get(ENV, "RADIX_ATDB", "")), "coheat.dat")
        isfile(coheat) || (coheat = joinpath(@__DIR__, "..", "data", "coheat.dat"))      # (the table of the repository is that of the XSTAR package)
        isfile(coheat) || (@info "Skipping the transfer tests (coheat.dat is needed)"; return)
        processes = Radix.standard_processes(Radix.load(Radix.Compton, coheat))
        mixture = Radix.Mixture(db, levels; multiplier=Dict(3 => 0.0, 4 => 0.0, 5 => 0.0))
        dir = joinpath(@__DIR__, "reference", "xstar_pow_xi2")
        spectrum = fits(joinpath(dir, "xout_cont1.fits"))[3].data
        lines = fits(joinpath(dir, "xout_lines1.fits"))[3].data
        edges = fits(joinpath(dir, "xout_rrc1.fits"))[3].data
        abundances = fits(joinpath(dir, "xout_abund1.fits"))[2].data
        E, L = Float64.(spectrum.energy), Float64.(spectrum.incident)*1e38
        run = Radix.iterate(mixture, 1e4, processes, E, L, [(1e13, 5e12), (1.5e13, 5e12)]; T=100.0, iterate=false)
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

        @slow "the continuum depth and the emission of a thick slab" begin
            # test/reference/xstar_thick_slab: 10²¹ cm⁻² at log ξ = 1, T = 10⁵ K, without the vturb of XSTAR (see its README). The zones are those of XSTAR: the rows
            # of its table are the states at their inner edges (`delta_r`), the last being the end of the slab
            dir = joinpath(@__DIR__, "reference", "xstar_thick_slab")
            slab = fits(joinpath(dir, "xout_cont1.fits"))[3].data
            Es, Is, Ts = Float64.(slab.energy), Float64.(slab.incident), Float64.(slab.transmitted)
            slabab = fits(joinpath(dir, "xout_abund1.fits"))[2].data
            r0 = 3.16228e20
            depth = Float64.(slabab.delta_r)[2:end]
            zones = [(r0 + depth[k], depth[k + 1] - depth[k]) for k in 1:length(depth) - 1]
            thickrun = Radix.iterate(mixture, 1e4, processes, Es, Is*1e38, zones; T=10.0, iterate=false, vturb=0.0, luminous=false)
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
            auto = Radix.slab(mixture, processes; ncn2=999, density=1e4, column=1e21, logξ=1.0, luminosity=1e8*1e38, α=-1.0, emult=1.0, taumax=5.0, steps=3, T=10.0, iterate=false, vturb=0.0, luminous=false)
            @test length(auto.zones) == 7 && auto.zones[1].Δr == 0
            @test auto.r ≈ r0 rtol=1e-6                              # the source of the reference run: E⁻¹ with the luminosity 10⁴⁶ erg/s
            @test auto.E ≈ Es rtol=1e-5
            @test auto.incident ./ 1e38 ≈ Is rtol=1e-5
            @test abs(sort(auto.transmitted ./ 1e38 ./ Ts)[end ÷ 2] - 1) < 5e-3 && all(abs.(auto.transmitted ./ 1e38 ./ Ts .- 1) .< 6e-2)
            @test [z.r - auto.r for z in auto.zones][1:6] ≈ [0; depth[1:5]] atol=3e-4*1.8e16 rtol=1e-3
            @test round.(log10.(1e4 .* cumsum([z.Δr for z in auto.zones])[2:end]); digits=2) == [20.26, 20.56, 20.74, 20.86, 20.95, 21.0]
            @test sum(z.Δr for z in auto.zones) ≈ 1e17
            # the final calculation of XSTAR, the last row of its table: a zone of no thickness at the end of the slab, with the temperature and electron fraction of the last zone
            @test (auto.final.r - auto.r) ≈ 1e17 rtol=1e-6
            @test auto.final.Δr == 0 && auto.final.T == auto.zones[end].T
            @test auto.zones[end].fractions[2][2] ≈ Float64(slabab.he_ii[end]) rtol=4e-3      # (the fractions of the last row of XSTAR's table are those of the last zone)
            @test auto.zones[end].fractions[findfirst(==(8), mixture.Z)][8] ≈ Float64(slabab.o_viii[end]) rtol=1e-3
            # three passes (`npass=3`: the outward optical depths of the lines and edges of the pass back from the outer edge): the ion fractions change by
            # 2.5% (He II in the first zone) to 0.2% and follow those of XSTAR in the change, which the escape probabilities of the outer side make
            auto3 = Radix.slab(mixture, processes; ncn2=999, density=1e4, column=1e21, logξ=1.0, luminosity=1e8*1e38, α=-1.0, emult=1.0, taumax=5.0, steps=3, T=10.0, iterate=false, vturb=0.0, luminous=false, passes=3)
            pass3 = fits(joinpath(dir, "xout_abund1_npass3.fits"))[2].data
            @test auto3.passes == 3 && length(auto3.zones) == 7
            for k in 1:7
                shift = auto3.zones[k].fractions[2][2]/auto.zones[k].fractions[2][2]
                @test shift ≈ Float64(pass3.he_ii[k])/Float64(slabab.he_ii[k]) atol=1.5e-3
                @test auto3.zones[k].fractions[2][2] ≈ Float64(pass3.he_ii[k]) rtol=4e-3
            end
            @test auto3.zones[1].fractions[2][2] < 0.98*auto.zones[1].fractions[2][2]
            @test_throws ArgumentError Radix.slab(mixture, processes; ncn2=999, density=1e4, column=1e21, logξ=1.0, luminosity=1e8*1e38, passes=2)
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
                @test all(r -> abs(r - 1) < 3.5e-2, bins)          # (5% before the recombination continua were emitted with the escape probabilities of their edges)
                @test sum(mine)/1e38 ≈ sum(reference) rtol=1e-2
            end
            # the end of the spectrum is the transmitted spectrum, where the lines do not absorb the bins that the edges do
            # the XSTAR of the reference does not absorb the continuum at all when vturb > 0 (gsmooth2)
            smoothed = fits(joinpath(dir, "xout_cont1_vturb1.fits"))[3].data
            @test count(==(1), Float64.(smoothed.transmitted) ./ Float64.(smoothed.incident)) > 700
        end

        @slow "the luminosities of the lines and the recombination edges of the thick slab" begin
            # the zones of `slab` (the seven of XSTAR) with its table of the lines (`xout_lines1.fits`: emit_inward + emit_outward in 10³⁸ erg/s, for the lines above 10²) and the total of
            # its table of the recombination edges (`xout_rrc1.fits`, 6 MB, not kept): 807150.6. The run of the table has the hydrogen of the database of the package of XSTAR
            slab = Radix.slab(mixture, processes; ncn2=999, density=1e4, column=1e21, logξ=1.0, luminosity=1e8*1e38, emult=1.0, taumax=5.0, steps=3, T=10.0, iterate=false, vturb=0.0)
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
            # the results as arrays and in one file
            model = Radix.Model(mixture, processes, Dict{Int, String}(Int(k) => v for (k, v) in ionlabel))
            results = Radix.tables(model, slab)
            @test results isa Radix.Results && propertynames(results) == (:zones, :spectrum, :lines, :edges)
            zone_table = results.zones
            @test length(zone_table.zone) == length(slab.zones) == 7
            @test zone_table.depth[1] == 0 && zone_table.thickness == [z.Δr for z in slab.zones] && issorted(zone_table.depth)
            @test zone_table.ion_parameter[1] ≈ 10 rtol=2e-5                                  # log ξ = 1 at the first zone
            @test zone_table.pressure ≈ zone_table.n_h .* 1.38e-12 .* zone_table.temperature rtol=1e-6
            @test zone_table.h_i .+ zone_table.h_ii ≈ ones(7)                                      # (the bare nucleus is named by the next stage)
            @test hasproperty(zone_table, :he_iii) && hasproperty(zone_table, :fe_xxvii) && hasproperty(zone_table, :zn_xxxi)
            @test all(z -> abs(sum(z.fractions[findfirst(==(2), mixture.Z)]) - 1) < 1e-9, slab.zones)
            spectrum_table = results.spectrum
            @test length(spectrum_table.energy) == 999 && spectrum_table.transmitted == slab.transmitted && spectrum_table.emit_inward == slab.spectrum.inward_continuum
            @test all(spectrum_table.transmitted_lines .<= spectrum_table.incident .* (1 + 1e-12))
            # the spectrum of xout_spect1: the transmitted spectrum with the lines and the emission (diffuse + the lines binned by `binemis`)
            spect = fits(joinpath(@__DIR__, "reference", "xstar_thick_slab", "xout_spect1.fits"))[3].data
            @test maximum(abs.(Float64.(spect.energy) ./ spectrum_table.energy .- 1)) < 1e-5
            visible = findall(>(1e-3*maximum(Float64.(spect.transmitted))), Float64.(spect.transmitted))              # (the bins that the file has more than 4 digits of)
            @test all(abs.(spectrum_table.transmitted_lines[visible] ./ 1e38 ./ Float64.(spect.transmitted)[visible] .- 1) .< 1e-3)
            for (diffuse, binned, reference) in ((spectrum_table.diffuse_inward, spectrum_table.lines_inward, spect.emit_inward), (spectrum_table.diffuse_outward, spectrum_table.lines_outward, spect.emit_outward))
                reference = Float64.(reference)
                strong = findall(>(1e-3*maximum(reference)), reference)
                ratio = (diffuse[strong] .+ binned[strong]) ./ 1e38 ./ reference[strong]
                @test length(strong) > 200
                @test abs(sort(ratio)[length(ratio) ÷ 2] - 1) < 2e-3
                @test count(r -> abs(r - 1) < 0.06, ratio) >= 0.93*length(ratio)          # (the others are hydrogen, the database of the package, and the bins of the continuum of He II)
                known !== nothing && known.name == "package" && @test count(r -> abs(r - 1) < 0.03, ratio) >= length(ratio) - 1          # (the package has the hydrogen of XSTAR: all the bins but one are within 3%; the one is 9.5% high, at 1.145 eV, where the He I triplet 10830 Å is: the populations of the He I of XSTAR are not converged, see the test of hydrogen and helium alone)
            end
            @test all(>=(0), spectrum_table.lines_inward) && all(>=(0), spectrum_table.lines_outward)
            @test sum(spectrum_table.lines_inward .* [Radix.xstar_energy_grid(999)[1]; diff(spectrum_table.energy)]) > 0
            line_table, edge_table = results.lines, results.edges
            @test length(line_table.ion) > 1000 && length(edge_table.ion) > 500
            @test sum(line_table.emit_inward) + sum(line_table.emit_outward) + sum(edge_table.emit_inward) + sum(edge_table.emit_outward) ≈ 1e38*(edge_total + sum(sum, slab.luminosities.inward[k][j] + slab.luminosities.outward[k][j]
                for k in eachindex(mixture.elements) for j in findall(c -> c isa Radix.AtomicLine2 && c.rtype == 4, mixture.elements[k].rates))) rtol=1e-3
            @test all(>(0), edge_table.energy)
            @test all(>(0), line_table.wavelength)
            @test "he_ii" in edge_table.ion && "he_ii" in line_table.ion
            path = tempname() * ".fits"
            write(path, results; meta=(; DENSITY=1e4, COLUMN=1e21))
            hdus = fits(path)
            @test length(hdus) == 5
            names = [strip(String(c.value)) for h in hdus[2:5] for c in h.cards if c.key == "EXTNAME"]
            @test names == ["ZONES", "SPECTRUM", "LINES", "EDGES"]
            back = hdus[2].data
            @test Float64.(back.temperature) ≈ zone_table.temperature rtol=1e-12
            @test Float64.(hdus[3].data.energy) == spectrum_table.energy
            @test strip.(hdus[4].data.ion) == line_table.ion && strip.(hdus[5].data.level) == edge_table.level
            @test count(c -> c.key == "TUNIT1", hdus[3].cards) == 1
            rm(path)
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

        @slow "a covering fraction of one half" begin
            # test/reference/xstar_cfrac: log ξ = 2 at 10⁴ cm⁻³ with the temperature iterated. The covering fraction cuts the photoexcitation of the lines, which XSTAR does in the balance and
            # the heating (the temperature of the first zone 4% lower: 18.081); the rows of its table are the first zone, again the first zone, and the states at the inner edges of the others
            if known !== nothing && known.name == "package"
                table = fits(joinpath(@__DIR__, "reference", "xstar_cfrac", "xout_abund1.fits"))[2].data
                r0 = Float64(table.radius[1])
                depth = Float64.(table.delta_r)[2:end]
                zones = [(r0 + depth[k], depth[k + 1] - depth[k]) for k in 1:length(depth) - 1]
                E = Radix.xstar_energy_grid(999)
                half = Radix.iterate(mixture, 1e4, Radix.standard_processes(Radix.load(Radix.Compton, coheat); cfrac=0.5), E, Radix.powerlaw(E, -1.0, 1e-2*1e38), zones; T=100.0, equilibrium=true,
                    vturb=0.0, luminous=false, cfrac=0.5)
                rows = [1, 3]
                @test [z.T for z in half.zones] ≈ Float64.(table.temperature)[rows] rtol=2e-3
                o = findfirst(==(8), mixture.Z)
                @test [z.fractions[o][8] for z in half.zones] ≈ Float64.(table.o_viii)[rows] rtol=5e-3
                @test [z.fractions[2][2] for z in half.zones] ≈ Float64.(table.he_ii)[rows] rtol=1e-2
            else
                @info "Skipping the test of the covering fraction (it is for the database of the XSTAR package)"
            end
        end

        @testset "the density from a file" begin
            # test/reference/xstar_density: XSTAR's `radexp=-100` (a line of r and n for each zone) at T = 10⁶ K. The zones are at the distances of the denseFile and the denseRun ends at the column or the end of the denseFile
            if known !== nothing && known.name == "package"
                denseDir = joinpath(@__DIR__, "reference", "xstar_density")
                denseFile = joinpath(denseDir, "density.dat")
                denseTable = Radix.DensityTable(denseFile)
                @test length(denseTable.r) == 40 && denseTable.n[1] == 7000.0 && Radix.gas_density(denseTable, 1.0) == 7000.0
                @test Radix.gas_density(denseTable, 1.0, denseTable.r[6]) == denseTable.n[6]
                @test Radix.gas_density(denseTable, 1.0, (denseTable.r[6] + denseTable.r[7])/2) ≈ (denseTable.n[6] + denseTable.n[7])/2
                @test_throws ArgumentError Radix.DensityTable([2.0, 1.0], [1.0, 1.0])
                denseO = findfirst(==(8), mixture.Z)
                denseKw = (; ncn2=999, luminosity=1e-3*1e38, α=-1.0, steps=2, emult=1.0, T=100.0, iterate=false, vturb=0.0, luminous=false)
                # the column of 10¹⁸ cm⁻²: five zones, and the final calculation at the sixth line
                denseRef = fits(joinpath(denseDir, "xout_abund1_column.fits"))[2].data
                denseRun = Radix.slab(mixture, processes; ncn2=999, density=denseFile, column=1e18, denseKw...)
                @test length(denseRun.zones) == 6 && denseRun.r == denseTable.r[1]
                @test [z.r - denseRun.r for z in denseRun.zones[1:5]] ≈ Float64.(denseRef.delta_r)[1:5] rtol=1e-5 atol=1e8
                @test (denseRun.final.r - denseRun.r) ≈ Float64(denseRef.delta_r[6]) rtol=1e-5
                @test denseRun.final.ntot == denseTable.n[7]
                @test [z.ntot for z in denseRun.zones[1:5]] ≈ Float64.(denseRef.n_p)[1:5] rtol=1e-5
                for k in 1:6          # (the last row has the fractions of the last zone)
                    @test denseRun.zones[k].fractions[2][2] ≈ Float64(denseRef.he_ii[k]) rtol=5e-3
                    @test denseRun.zones[k].fractions[denseO][8] ≈ Float64(denseRef.o_viii[k]) rtol=5e-3
                end
                # the final calculation (T and x_e of the last zone, its radiation, the density and distance of the line): the heating and cooling of the photons in XSTAR's denseTable
                @test denseRun.final.heating ≈ 1.14746613e-15 rtol=3e-3
                @test denseRun.final.cooling ≈ 1.26292604e-15 rtol=3e-3
                # the column of 10¹⁹ cm⁻² is more than the file has: the zones are its lines, and the end of the slab is the last
                denseRef = fits(joinpath(denseDir, "xout_abund1_table.fits"))[2].data
                denseRun = Radix.slab(mixture, processes; ncn2=999, density=denseFile, column=1e19, denseKw...)
                @test length(denseRun.zones) == 40 && denseRun.final.r == denseTable.r[end]
                @test [z.r - denseRun.r for z in denseRun.zones] ≈ Float64.(denseRef.delta_r) rtol=1e-5 atol=1e8
                @test [z.fractions[2][2] for z in denseRun.zones] ≈ Float64.(denseRef.he_ii) rtol=5e-3
                @test [z.fractions[denseO][8] for z in denseRun.zones] ≈ Float64.(denseRef.o_viii) rtol=5e-3
            else
                @info "Skipping the test of the density from a file (it is for the database of the XSTAR package)"
            end
        end

        @testset "the thermal sources with the temperature iterated" begin
            # test/reference/xstar_sources: a blackbody (kT = 431 eV) and thermal bremsstrahlung (kT = 86 eV) at log ξ = 2, 10⁴ cm⁻³, 10³⁵ erg/s, T and x_e iterated (a column of
            # 10¹⁷ cm⁻²: the first zone and the end of the slab; the rows of the table are the first zone, again the first zone, and the end). For the database of the package only
            if known !== nothing && known.name == "package"
                oxygen = findfirst(==(8), mixture.Z)
                for (name, spectrum, tolerance) in (("bbody", (E, L) -> Radix.blackbody(E, 0.5, L), 2e-3), ("brems", (E, L) -> Radix.thermal(E, 0.1, L), 6e-3))
                    sourcetable = fits(joinpath(@__DIR__, "reference", "xstar_sources", "thin_$(name)_abund1.fits"))[2].data
                    sourcerun = Radix.slab(mixture, processes; ncn2=999, density=1e4, column=1e17, logξ=2.0, luminosity=1e-3*1e38, spectrum, steps=2, emult=0.5, T=100.0, equilibrium=true, vturb=0.0, luminous=false)
                    sourcezone = sourcerun.zones[1]
                    @test sourcezone.T ≈ Float64(sourcetable.temperature[1]) rtol=tolerance
                    @test sourcezone.xee ≈ Float64(sourcetable.x_e[1]) rtol=1e-4
                    @test sourcezone.fractions[2][2] ≈ Float64(sourcetable.he_ii[1]) rtol=4e-3
                    @test sourcezone.fractions[oxygen][8] ≈ Float64(sourcetable.o_viii[1]) rtol=4e-3
                    @test sourcezone.fractions[1][1] ≈ Float64(sourcetable.h_i[1]) rtol=4e-3
                end
            else
                @info "Skipping the test of the thermal sources (it is for the database of the XSTAR package)"
            end
        end

        @slow "the thermal sources through a slab" begin
            # the same sources through 10¹⁹ cm⁻² (emult = 1: the zones of the table are those of Radix, which agree to 10⁻⁴ in r). The first four rows (T from 11.8 to 7.6 in 10⁴ K for the
            # blackbody) agree to 4×10⁻³ in T and 1% in the fractions, and the brems slab (to 3.4) throughout. The blackbody slab ends in a thermal front, where T falls to 1.8 within
            # one zone and the trace ions of C to Fe change by factors (XSTAR's own are not converged there: its T depends by 1% on the steps): T agrees to 5×10⁻³ and 2×10⁻² at the end
            if known !== nothing && known.name == "package"
                oxygen = findfirst(==(8), mixture.Z)
                for (name, spectrum, rows, tolerance, front) in (("bbody", (E, L) -> Radix.blackbody(E, 0.5, L), 1:4, 4e-3, true), ("brems", (E, L) -> Radix.thermal(E, 0.1, L), 1:5, 7e-3, false))
                    sourcetable = fits(joinpath(@__DIR__, "reference", "xstar_sources", "thick_$(name)_abund1.fits"))[2].data
                    sourcerun = Radix.slab(mixture, processes; ncn2=999, density=1e4, column=1e19, logξ=2.0, luminosity=1e-3*1e38, spectrum, steps=2, emult=1.0, T=100.0, equilibrium=true, vturb=0.0, luminous=false)
                    sourcestates = [sourcerun.zones[1:end - 1]..., sourcerun.final]            # (XSTAR's table has no row for the last zone, the final calculation is its last)
                    @test length(sourcestates) == length(sourcetable.temperature)
                    @test [z.r - sourcerun.r for z in sourcestates] ≈ Float64.(sourcetable.delta_r) rtol=2e-4 atol=1e10
                    for k in rows
                        @test sourcestates[k].T ≈ Float64(sourcetable.temperature[k]) rtol=tolerance
                        @test sourcestates[k].fractions[2][2] ≈ Float64(sourcetable.he_ii[k]) rtol=1e-2
                        @test sourcestates[k].fractions[oxygen][8] ≈ Float64(sourcetable.o_viii[k]) rtol=1e-2
                        @test sourcestates[k].fractions[1][1] ≈ Float64(sourcetable.h_i[k]) rtol=1.5e-2
                    end
                    if front
                        @test sourcestates[5].T ≈ Float64(sourcetable.temperature[5]) rtol=5e-3
                        @test sourcestates[end].T ≈ Float64(sourcetable.temperature[end]) rtol=2e-2
                        @test sourcerun.zones[end].fractions[2][2] ≈ Float64(sourcetable.he_ii[end]) rtol=1e-2      # (the fractions of the last row are the last zone's)
                    end
                end
            else
                @info "Skipping the test of the thermal sources through a slab (it is for the database of the XSTAR package)"
            end
        end

        @testset "hydrogen and helium alone: the populations of the He I of XSTAR are not converged" begin
            # test/reference/xstar_hhe: the thick slab with the other elements left out, whose table of the lines lists them all. The lines of H I, He II and the singlets of He I agree; the triplets of He I
            # (10830 Å, ...) are 2 times weaker in XSTAR because its populations of the He I are not a steady state of its own matrix (the relaxation of `msolvelucy` leaves out the levels below 10⁻⁶)
            if known !== nothing && known.name == "package"
                hhe_dir = joinpath(@__DIR__, "reference", "xstar_hhe")
                hhemix = Radix.Mixture(db, levels; multiplier=Dict(Z => 0.0 for Z in 3:30))
                hhemodel = Radix.Model(hhemix, processes, Dict{Int, String}(Int(r.ion) => strip(r.label) for r in db if r isa Radix.Ion))
                hherun = Radix.slab(hhemix, processes; ncn2=999, density=1e4, column=1e21, logξ=1.0, luminosity=1e8*1e38, emult=1.0, taumax=5.0, steps=3, T=10.0, iterate=false, vturb=0.0)
                hhelines = Radix.tables(hhemodel, hherun).lines
                hhe_mine = Dict{Tuple{String, String, String}, Float64}()
                for i in eachindex(hhelines.ion)
                    hhe_key = (strip(hhelines.ion[i]), strip(hhelines.lower[i]), strip(hhelines.upper[i]))
                    hhe_mine[hhe_key] = get(hhe_mine, hhe_key, 0.0) + (hhelines.emit_inward[i] + hhelines.emit_outward[i])/1e38
                end
                hhe_rows = [split(line, ',') for line in readlines(joinpath(hhe_dir, "lines.csv"))[2:end]]
                triplet(hhe_row) = hhe_row[1] == "he_i" && occursin(r"\.3[SPDFG]", hhe_row[2] * hhe_row[3])
                hhe_singles, hhe_triples = Float64[], Float64[]
                for hhe_row in hhe_rows
                    hhe_emission = parse(Float64, hhe_row[5]) + parse(Float64, hhe_row[6])
                    hhe_emission > 1 || continue
                    hhe_key = (hhe_row[1], hhe_row[2], hhe_row[3])
                    hhe_ratio = get(hhe_mine, hhe_key, get(hhe_mine, (hhe_row[1], hhe_row[3], hhe_row[2]), NaN))/hhe_emission
                    push!(triplet(hhe_row) ? hhe_triples : hhe_singles, hhe_ratio)
                end
                @test length(hhe_singles) > 150 && count(r -> abs(r - 1) < 0.05, hhe_singles) >= 0.97*length(hhe_singles) && abs(sort(hhe_singles)[end ÷ 2] - 1) < 5e-3
                @test length(hhe_triples) >= 3 && all(r -> 1.7 < r < 2.4, hhe_triples)
                # XSTAR's populations in the matrix of the first zone: He I 2³S gets 1.5 times what it loses, and Radix's solution is the steady state
                hhe_zone = hherun.zones[1]
                hhe_k = findfirst(==(2), hhemix.Z)
                hhe_A = Matrix(Radix.element_matrix(hhemix.elements[hhe_k], Radix.Cell(10.0, hhe_zone.nₕ, hhe_zone.nₑ, hhe_zone.ntot); radiation=hhe_zone.radiation, escape=hhe_zone.escape[hhe_k]))
                hhe_pops = [split(line, ',') for line in readlines(joinpath(hhe_dir, "populations_zone1.csv"))[2:end]]
                hhe_hei = [parse(Float64, p[3]) for p in hhe_pops if p[1] == "he_i"]
                hhe_heii = [parse(Float64, p[3]) for p in hhe_pops if p[1] == "he_ii"]
                hhe_x = vcat(hhe_hei[1:45], hhe_heii)                    # (the continuum of He I is the ground state of He II)
                hhe_inflow = sum(hhe_A[2, j]*hhe_x[j] for j in eachindex(hhe_x) if j != 2)
                @test hhe_inflow > 1.4*(-hhe_A[2, 2]*hhe_x[2])
                @test hhe_zone.populations[hhe_k][2] > 4*hhe_x[2]
                @test hhe_zone.populations[hhe_k][2] ≈ Radix.level_populations(hhe_A)[2] rtol=1e-6
            else
                @info "Skipping the test of hydrogen and helium alone (it is for the database of the XSTAR package)"
            end
        end

        @slow "three passes: the spectrum, the lines and the edges" begin
            # test/reference/xstar_thick_slab_npass3: the thick slab with `npass=3`. The continuum, the binned emission, the luminosities of the lines (and how they change from the first pass: the He II 3s-2p line by 1.6)
            # and the total of the recombination edges follow XSTAR's (package database only)
            if known !== nothing && known.name == "package"
                p3dir = joinpath(@__DIR__, "reference", "xstar_thick_slab_npass3")
                p3model = Radix.Model(mixture, processes, Dict{Int, String}(Int(r.ion) => strip(r.label) for r in db if r isa Radix.Ion))
                p3kw = (; ncn2=999, density=1e4, column=1e21, logξ=1.0, luminosity=1e8*1e38, emult=1.0, taumax=5.0, steps=3, T=10.0, iterate=false, vturb=0.0)
                p3one = Radix.tables(p3model, Radix.slab(mixture, processes; ncn2=999, p3kw..., passes=1))
                p3run = Radix.slab(mixture, processes; ncn2=999, p3kw..., passes=3)
                p3res = Radix.tables(p3model, p3run)
                p3sp = p3res.spectrum
                p3cont = fits(joinpath(p3dir, "xout_cont1.fits"))[3].data
                p3spect = fits(joinpath(p3dir, "xout_spect1.fits"))[3].data
                p3t = Float64.(p3cont.transmitted)
                p3vis = findall(>(1e-3*maximum(p3t)), p3t)
                @test all(abs.((p3sp.transmitted[p3vis] ./ 1e38) ./ p3t[p3vis] .- 1) .< 1e-3)
                for (p3mine, p3ref) in ((p3sp.emit_inward, p3cont.emit_inward), (p3sp.emit_outward, p3cont.emit_outward))
                    p3ref = Float64.(p3ref); p3strong = findall(>(1e-3*maximum(p3ref)), p3ref)
                    p3b = (p3mine[p3strong] ./ 1e38) ./ p3ref[p3strong]
                    @test abs(sort(p3b)[length(p3b) ÷ 2] - 1) < 2e-3 && all(b -> abs(b - 1) < 3.5e-2, p3b)
                end
                for (p3d, p3l, p3ref) in ((p3sp.diffuse_inward, p3sp.lines_inward, p3spect.emit_inward), (p3sp.diffuse_outward, p3sp.lines_outward, p3spect.emit_outward))
                    p3ref = Float64.(p3ref); p3strong = findall(>(1e-3*maximum(p3ref)), p3ref)
                    p3b = (p3d[p3strong] .+ p3l[p3strong]) ./ 1e38 ./ p3ref[p3strong]
                    @test length(p3strong) > 200 && abs(sort(p3b)[length(p3b) ÷ 2] - 1) < 2e-3
                    @test count(b -> abs(b - 1) < 0.07, p3b) >= length(p3b) - 3             # (the He I bin at 1.145 eV, the C V forbidden line at 41 Å and a bin at 245 Å)
                end
                # the luminosities of the lines, and their change from the first pass
                p3agg(t) = (d = Dict{Tuple{String, String, String}, Float64}(); for i in eachindex(t.ion); p3k = (strip(t.ion[i]), strip(t.lower[i]), strip(t.upper[i])); d[p3k] = get(d, p3k, 0.0) + (t.emit_inward[i] + t.emit_outward[i])/1e38; end; d)
                p3m3, p3m1 = p3agg(p3res.lines), p3agg(p3one.lines)
                p3x3 = fits(joinpath(p3dir, "xout_lines1.fits"))[3].data
                p3x1 = fits(joinpath(p3dir, "xout_lines1_npass1.fits"))[3].data
                p3x1d = Dict{Tuple{String, String, String}, Float64}(); for i in eachindex(p3x1.ion); p3x1d[(strip(p3x1.ion[i]), strip(p3x1.lower_level[i]), strip(p3x1.upper_level[i]))] = Float64(p3x1.emit_inward[i]) + Float64(p3x1.emit_outward[i]); end
                p3ratio, p3change = Float64[], Float64[]
                for i in eachindex(p3x3.ion)
                    p3k = (strip(p3x3.ion[i]), strip(p3x3.lower_level[i]), strip(p3x3.upper_level[i]))
                    p3v = Float64(p3x3.emit_inward[i]) + Float64(p3x3.emit_outward[i])
                    p3a3 = get(p3m3, p3k, get(p3m3, (p3k[1], p3k[3], p3k[2]), NaN)); p3a1 = get(p3m1, p3k, get(p3m1, (p3k[1], p3k[3], p3k[2]), NaN))
                    isfinite(p3a3) && push!(p3ratio, p3a3/p3v)
                    (haskey(p3x1d, p3k) && isfinite(p3a3) && isfinite(p3a1)) && push!(p3change, (p3a3/p3a1)/(p3v/p3x1d[p3k]))
                end
                @test length(p3ratio) > 590 && abs(sort(p3ratio)[end ÷ 2] - 1) < 3e-3 && count(r -> abs(r - 1) < 0.05, p3ratio) >= 0.98*length(p3ratio)
                @test length(p3change) > 540 && count(c -> abs(c - 1) < 0.02, p3change) >= 0.95*length(p3change)           # (XSTAR changes 414 of them by more than 1%)
                p3he = ("he_ii", "1s0.2p1.2P_1/2", "1s0.3s1.2S_1/2")
                @test p3m3[p3he]/p3m1[p3he] ≈ 1.596 rtol=0.03
                # the total of the recombination edges (XSTAR's table: 798877.05, and 807150.61 with one pass)
                p3e = p3res.edges
                @test sum(p3e.emit_inward .+ p3e.emit_outward)/1e38 ≈ 798877.05 rtol=3e-3
            else
                @info "Skipping the test of the three passes (it is for the database of the XSTAR package)"
            end
        end

        @testset "levels below 1e-6: the two lines that XSTAR changes between passes" begin
            # test/reference/xstar_hcar: hydrogen, helium, carbon and argon only. The forbidden 2³S line of C V and the 3s ¹P line of Ar IX change by +10.5% and +13.7% from the first pass of XSTAR to the third (Radix: 0 and +4%),
            # because the populations of their upper levels (below 10⁻⁶) are not a steady state of XSTAR's own matrix: XSTAR's populations in the first and the third zones against the matrix of Radix
            if known !== nothing && known.name == "package"
                hcar_dir = joinpath(@__DIR__, "reference", "xstar_hcar")
                hcar_mix = Radix.Mixture(db, levels; multiplier=Dict(Z => 0.0 for Z in 3:30 if Z != 6 && Z != 18))
                hcar_label = Dict{Int, String}(Int(r.ion) => lowercase(strip(r.label)) for r in db if r isa Radix.Ion)
                hcar_model = Radix.Model(hcar_mix, processes, Dict{Int, String}(Int(r.ion) => strip(r.label) for r in db if r isa Radix.Ion))
                hcar_kw = (; ncn2=999, density=1e4, column=1e21, logξ=1.0, luminosity=1e8*1e38, emult=1.0, taumax=5.0, steps=3, T=10.0, iterate=false, vturb=0.0)
                hcar_one = Radix.slab(hcar_mix, processes; ncn2=999, hcar_kw..., passes=1)
                hcar_three = Radix.slab(hcar_mix, processes; ncn2=999, hcar_kw..., passes=3)
                hcar_lines(run) = (t = Radix.tables(hcar_model, run).lines; d = Dict{Tuple{String, String, String}, Float64}(); for i in eachindex(t.ion); d[(strip(t.ion[i]), strip(t.lower[i]), strip(t.upper[i]))] = (t.emit_inward[i] + t.emit_outward[i])/1e38; end; d)
                hcar_x(file) = (d = Dict{Tuple{String, String, String}, Float64}(); for line in readlines(joinpath(hcar_dir, file))[2:end]; f = split(line, ','); d[(f[1], f[2], f[3])] = parse(Float64, f[5]) + parse(Float64, f[6]); end; d)
                hcar_r1, hcar_r3, hcar_x1, hcar_x3 = hcar_lines(hcar_one), hcar_lines(hcar_three), hcar_x("lines_npass1.csv"), hcar_x("lines_npass3.csv")
                hcar_cv = ("c_v", "1s2.1S_0", "1s1.2s1.3S_1"); hcar_ar = ("ar_ix", "2p6.1S_0", "2p5.3s1.1P_1")
                # the C V triplet 2³P-3³D at 248.7 Å drops by 16% in XSTAR from the first pass to the third (unconverged populations; Radix: -0.06%) below a line of 257 that Radix
                # ranks 11th in its bin, which changes which 10 lines of the bin are spread (the 245 Å bin of the three passes of the thick slab)
                hcar_tr = ("c_v", "1s1.2p1.3P_1", "1s1.3d1.3D")
                @test hcar_x3[hcar_tr]/hcar_x1[hcar_tr] ≈ 0.843 rtol=5e-3
                @test abs(hcar_r3[hcar_tr]/hcar_r1[hcar_tr] - 1) < 2e-3
                @test hcar_x3[hcar_cv]/hcar_x1[hcar_cv] ≈ 1.105 rtol=5e-3
                @test abs(hcar_r3[hcar_cv]/hcar_r1[hcar_cv] - 1) < 5e-3
                @test hcar_x3[hcar_ar]/hcar_x1[hcar_ar] ≈ 1.137 rtol=5e-3
                @test 1.02 < hcar_r3[hcar_ar]/hcar_r1[hcar_ar] < 1.07
                hcar_resid(Z, zone, file, level) = begin
                    k = findfirst(==(Z), hcar_mix.Z); el = hcar_mix.elements[k]; z = hcar_one.zones[zone]
                    A = Matrix(Radix.element_matrix(el, Radix.Cell(10.0, z.nₕ, z.nₑ, z.ntot); radiation=z.radiation, escape=z.escape[k]))
                    pop = Dict{Tuple{String, String}, Float64}()
                    for line in readlines(joinpath(hcar_dir, file))[2:end]; f = split(line, ','); haskey(pop, (f[1], f[2])) || (pop[(f[1], f[2])] = parse(Float64, f[3])); end
                    x = zeros(size(A, 1)); g = 0
                    for (j, ion) in enumerate(el.ions), l in 1:(j == length(el.ions) ? el.nlev[j] : el.nlev[j] - 1)
                        key = (hcar_label[ion], strip(el.levels[(ion, l)].label))
                        haskey(pop, key) && (x[el.offset[j] + l] = pop[key])
                        key == level && (g = el.offset[j] + l)
                    end
                    inflow = sum(max(A[g, j]*x[j], 0.0) for j in eachindex(x) if j != g); outflow = -A[g, g]*x[g]
                    ((inflow - outflow)/max(inflow, outflow), z.populations[k][g]/x[g])
                end
                hcar_a = hcar_resid(6, 1, "populations_zone1.csv", ("c_v", "1s1.2s1.3S_1"))
                @test hcar_a[1] < -0.9 && hcar_a[2] < 0.2                       # (the 2³S of C V in the first zone: 6% of what it loses is what it gets, and Radix's population is 7 times lower)
                hcar_b = hcar_resid(6, 3, "populations_zone3.csv", ("c_v", "1s1.2s1.3S_1"))
                @test -0.06 < hcar_b[1] < -0.01 && 0.7 < hcar_b[2] < 0.9
                hcar_c = hcar_resid(18, 3, "populations_zone3.csv", ("ar_ix", "2p5.3s1.3P_2"))
                @test hcar_c[1] < -0.05 && 0.7 < hcar_c[2] < 0.85                   # (and the first zone of the Ar IX agrees to 0.2%)
                hcar_d = hcar_resid(18, 1, "populations_zone1.csv", ("ar_ix", "2p5.3s1.3P_2"))
                @test abs(hcar_d[1]) < 2e-3 && abs(hcar_d[2] - 1) < 5e-3
            else
                @info "Skipping the test of the levels below 1e-6 (it is for the database of the XSTAR package)"
            end
        end

        @slow "10000 energy bins: the balance is on the 999 bins of XSTAR, the continuum arrays on the 10000" begin
            # test/reference/xstar_10000bins: XSTAR's balance, heating and edges (`dsec`, `calc_hmc_all`, `calc_emisab_all`) use the field that `bremsmap` maps to 999 bins, whatever `ncn2` is: the temperature of the thin slab
            # changes by 0.03% from 999 to 10000 bins in XSTAR (18.8925 and 18.8982) and Radix follows (18.8958 and 18.9015, against 18.71 when the balance was on the bins of the field). The edges that `rlbin` ranks
            # include those with an emissivity but no opacity, which the bins of 10000 (fewer edges in each) let in, so the continuum emission of the thick slab agrees with that of XSTAR
            if known !== nothing && known.name == "package"
                n10dir = joinpath(@__DIR__, "reference", "xstar_10000bins")
                n10thin = fits(joinpath(n10dir, "thin_abund1.fits"))[2].data
                n10o = findfirst(==(8), mixture.Z)
                for (n10n, n10tol) in ((999, 2e-3), (10000, 2e-3))
                    n10run = Radix.slab(mixture, processes; density=1e4, column=1e17, logξ=2.0, luminosity=1e-6*1e38, T=100.0, equilibrium=true, vturb=0.0, luminous=false, steps=2, emult=0.5, ncn2=n10n)
                    @test n10run.zones[1].T ≈ Float64(n10thin.temperature[1]) rtol=n10tol
                    n10n == 10000 && @test n10run.zones[1].fractions[n10o][8] ≈ Float64(n10thin.o_viii[1]) rtol=5e-3
                end
                n10slab = fits(joinpath(n10dir, "thick_abund1.fits"))[2].data
                n10cont = fits(joinpath(n10dir, "thick_cont1.fits"))[3].data
                n10E, n10I = Float64.(n10cont.energy), Float64.(n10cont.incident)
                n10r0 = 3.16228e20
                n10depth = Float64.(n10slab.delta_r)[2:end]
                n10zones = [(n10r0 + n10depth[k], n10depth[k + 1] - n10depth[k]) for k in 1:length(n10depth) - 1]
                n10thick = Radix.iterate(mixture, 1e4, processes, n10E, n10I*1e38, n10zones; T=10.0, iterate=false, vturb=0.0, luminous=false)
                n10t = Float64.(n10cont.transmitted); n10vis = findall(>(1e-3*maximum(n10t)), n10t)
                @test all(abs.(exp.(-n10thick.dpthcont[n10vis]) .* n10I[n10vis] ./ n10t[n10vis] .- 1) .< 1e-3)
                for (n10mine, n10ref) in ((n10thick.spectrum.inward_continuum, n10cont.emit_inward), (n10thick.spectrum.outward_continuum, n10cont.emit_outward))
                    n10ref = Float64.(n10ref); n10strong = findall(>(1e-3*maximum(n10ref)), n10ref)
                    n10b = (n10mine[n10strong] ./ 1e38) ./ n10ref[n10strong]
                    @test length(n10strong) > 5000 && abs(sort(n10b)[length(n10b) ÷ 2] - 1) < 3e-3
                    @test count(b -> abs(b - 1) > 0.03, n10b) <= 12 && all(b -> abs(b - 1) < 0.06, n10b)
                    @test sum(n10mine)/1e38/sum(n10ref) ≈ 1 rtol=4e-3
                end
                for (n10k, n10row) in ((1, 2), (6, 7), (13, 14))
                    @test n10thick.zones[n10k].fractions[2][2] ≈ Float64(n10slab.he_ii[n10row]) rtol=1e-2
                end
                # XSTAR keeps the ten strongest edges and lines of each bin (`rank_per_bin`, 10 in the constants of `ucalc()` that these tests use): the continuum emission then depends on the number of bins
                # (13% in the total from 999 to 30000 bins), and with the 1000 of the default constants it hardly does
                n10total(n) = begin
                    n10grid = Radix.xstar_energy_grid(n)
                    n10run = Radix.iterate(mixture, 1e4, processes, n10grid, Radix.powerlaw(n10grid, -1.0, 1e8*1e38), n10zones; T=10.0, iterate=false, vturb=0.0, luminous=false)
                    n10f = n10run.spectrum.inward_continuum
                    n10w(f, lo, hi) = sum(((1 - f[i]) + (1 - f[i + 1]))/2*(n10grid[i + 1] - n10grid[i]) for i in 1:length(n10grid) - 1 if lo <= n10grid[i] < hi)
                    (sum((n10f[i] + n10f[i + 1])/2*(n10grid[i + 1] - n10grid[i]) for i in 1:length(n10grid) - 1), n10w(n10run.transmission, 13.6, 54.0), n10w(exp.(-n10run.dpthc), 13.6, 54.0))
                end
                n10xstar = n10total(999), n10total(10000)
                n10converged = Radix.constants(Radix.Constants(Radix.constants(); rank_per_bin=Radix.converged_rank_per_bin, exact_binning=true, resolved_lines=true)) do
                    n10total(999), n10total(10000)
                end
                @test abs(n10xstar[2][1]/n10xstar[1][1] - 1) > 0.04
                @test abs(n10converged[2][1]/n10converged[1][1] - 1) < 0.02
                # the equivalent width of the absorption of the slab in 14-54 eV with its lines: the mean opacity of the bins (exp(-dpthc), which XSTAR's transmitted spectrum with lines is) gives 2 to 3 times as much at 999 bins
                # as at 10000, and the profiles of the lines inside the bins (`resolved_transmission`) 7%
                @test n10converged[1][3]/n10converged[2][3] > 1.5
                @test abs(n10converged[1][2]/n10converged[2][2] - 1) < 0.12
            else
                @info "Skipping the test of the 10000 energy bins (it is for the database of the XSTAR package)"
            end
        end

        @slow "three passes with the temperature iterated" begin
            # test/reference/xstar_thick_slab/xout_abund1_T_npass{1,3}.fits: the thick slab of the run above with niter=99 (T and x_e iterated), npass=1 and 3. The zones are those of
            # XSTAR (its rows are the states at the inner edges: the first zone is row 1 and the others the row after their number). XSTAR's own temperatures are those of
            # the database of its package (the one of the tree has other hydrogen records: 3-6% in T), so this is for that database only. The three passes raise the temperature of
            # the zones by 0.9% (4.5513 to 4.5927 in the first), and XSTAR's `dsec` leaves its first rows of one pass at 4.5513 (not converged), which Radix does not
            if known !== nothing && known.name == "package"
                dir = joinpath(@__DIR__, "reference", "xstar_thick_slab")
                table = fits(joinpath(dir, "xout_abund1_T_npass3.fits"))[2].data
                depth = Float64.(table.delta_r)[2:end]
                r0 = 3.16228e20
                zones = [(r0 + depth[k], depth[k + 1] - depth[k]) for k in 1:length(depth) - 1]
                E = Radix.xstar_energy_grid(999)
                three = Radix.iterate(mixture, 1e4, processes, E, Radix.powerlaw(E, -1.0, 1e8*1e38), zones; T=10.0, equilibrium=true, vturb=0.0, luminous=false, passes=3)
                rows = [1; 3:8]
                T = Float64.(table.temperature)[rows]
                @test [z.T for z in three.zones[1:7]] ≈ T rtol=2e-3
                he_ii, o_viii = Float64.(table.he_ii)[rows], Float64.(table.o_viii)[rows]
                @test [z.fractions[2][2] for z in three.zones[1:7]] ≈ he_ii rtol=1.5e-2
                o = findfirst(==(8), mixture.Z)
                @test [z.fractions[o][8] for z in three.zones[1:7]] ≈ o_viii rtol=5e-3
                one = Float64.(fits(joinpath(dir, "xout_abund1_T_npass1.fits"))[2].data.temperature)
                @test all(three.zones[k].T > one[k + 1] for k in 4:7)                # (the passes raise the temperature)
            else
                @info "Skipping the test of three passes with the temperature iterated (it is for the database of the XSTAR package)"
            end
        end

        @testset "a slab at constant pressure" begin
            # test/reference/xstar_const_pressure: `lcpres=1 pressure=1e-7`, log Ξ = 1, 10²¹ cm⁻² (the thick slab at a pressure instead of a density): r = 1.62965×10²⁰ cm, the first zone at
            # T = 10.4611 and n = 6927 (the database of the package of XSTAR: 10.4574 and 6929 here, 0.04% and 0.03%; the one of the tree has other hydrogen records)
            pressure = 1e-7
            cp = Radix.slab(mixture, processes; ncn2=999, pressure, column=3e20, logξ=1.0, luminosity=1e8*1e38, emult=1.0, taumax=5.0, steps=3, T=10.0, equilibrium=true, vturb=0.0, luminous=false)
            @test cp.r ≈ 1.62965e20 rtol=1e-5
            @test length(cp.zones) == 2 && cp.zones[1].Δr == 0
            @test all(z -> isapprox(z.ntot*1.38e-12*z.T, pressure; rtol=1e-6), cp.zones)          # n = P/(k T) at the temperature of each zone
            @test cp.zones[1].T > 8 && cp.zones[1].T < 12 && 6e3 < cp.zones[1].ntot < 8e3
            @test sum(z.ntot*z.Δr for z in cp.zones) ≈ 3e20 rtol=1e-9                           # the column of the zones
            @test all(z -> z.converged && abs(z.imbalance) <= 1e-4, cp.zones)
        end

        @testset "a density that falls as a power of the distance" begin
            # `radexp=-3` (a run of XSTAR with nsteps=3, emult=1, rlrad38=1e-2, log ξ = 1, T = 10⁵ K kept: the column of n ∝ r⁻³ from 3.16×10¹⁵ cm is limited, so that the run never ends; its log shows the rows
            # log r = 15.50, 15.50, 15.62, 15.75, 15.87, 16.00; log N = -10, -10, 18.65, 18.84, 18.92, 18.96; log n = 4.00, 4.00, 3.63, 3.25, 2.88, 2.50): the zones of r/3 of its first rows
            r0 = Radix.source_distance(1e-2*1e38, 10.0, 1e4)
            @test log10(r0) ≈ 15.5 atol=0.005
            density = Radix.PowerLawDensity(1e4, r0, -3.0)
            zones = [(r0, 0.0), (r0, r0/3), (4r0/3, 4r0/9)]
            rad = Radix.iterate(mixture, density, processes, Radix.xstar_energy_grid(999), Radix.powerlaw(Radix.xstar_energy_grid(999), -1.0, 1e-2*1e38), zones; T=10.0,
                iterate=false, vturb=0.0, luminous=false)
            @test round.(log10.([z.r for z in rad.zones]); digits=2) == [15.5, 15.5, 15.62]
            @test round.(log10.([z.ntot for z in rad.zones]); digits=2) == [4.0, 4.0, 3.63]       # (the density of the inner edge of each zone)
            column = cumsum([Radix.gas_density(density, 10.0, z.r + z.Δr)*z.Δr for z in rad.zones])        # (added with the density of the outer edge, as XSTAR does)
            @test round(log10(column[2]); digits=2) == 18.65
            @test round(log10(column[3]); digits=2) == 18.84
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

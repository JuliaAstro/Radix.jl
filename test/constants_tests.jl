# The physical constants of the rates: CODATA 2022 from PhysicalConstants.jl by default, the rounded
# values of XSTAR's ucalc on request.

@testset "Physical constants" begin
    K = Radix.constants()
    U = Radix.ucalc_constants()
    CODATA2018 = Radix.PhysicalConstants.CODATA2018

    @testset "CODATA 2022 is the default" begin
        @test K == Radix.Constants()
        @test K.ergsev == 1.602176634e-12             # exact in the SI
        @test K.kB_cgs == 1.380649e-16
        @test K.Ry_eV ≈ 13.605693122990 rtol=1e-12
        @test K.hc_eVÅ ≈ 12398.419843320 rtol=1e-12
        @test K.hc_over_k ≈ 1.438776877e8 rtol=1e-9
        @test K.proton_mass ≈ 1836.152673 rtol=1e-8
        # π is exact; ucalc's 4π, 8π and π are rounded
        @test K.fourpi == Float64(4π) && K.eightpi == Float64(8π) && K.pi_pexs == Float64(π)
        @test U.fourpi == 12.56 && U.eightpi == 25.3 && U.pi_pexs == 3.14159
    end

    @testset "consistency of the derived constants" begin
        @test K.T_floor_coeff == K.hc_over_k/Radix.T_floor_dE_over_kT
        @test K.hc_over_k ≈ K.hc_eVÅ/K.kB_eV rtol=1e-12
        @test K.kT_eV == K.kB_eV*Radix.T_unit
        @test K.Ry_K ≈ K.Ry_eV/K.kB_eV rtol=1e-12
        @test K.eV_K ≈ 1/K.kB_eV rtol=1e-12
        @test K.inv_kB*K.kB_cgs ≈ 1 rtol=1e-12
        @test K.Ry_erg ≈ K.Ry_eV*K.ergsev rtol=1e-12
        @test K.invcm_per_eV*K.hc_eVÅ ≈ 1e8 rtol=1e-12      # cm⁻¹ eV⁻¹ × eV Å = Å per cm
        @test K.Rinf_invcm ≈ K.Ry_eV*K.invcm_per_eV rtol=1e-12
        # the one physical constant ucalc spells several ways is one value here
        @test K.ergsev_bremsint == K.ergsev == K.ergsev_decay
        @test K.kB_szirc == K.kB_cgs
        @test K.Ry_K == K.Ry_K_erc == K.Ry_K_cb
        @test K.Ry_eV_single == K.Ry_eV && K.hc_eVÅ_single == K.hc_eVÅ
    end

    @testset "ucalc's rounded values are close" begin
        # every field within the rounding of XSTAR's constants; the known larger deviations are the
        # proton mass (ucalc takes 1800 electron masses), 8π (25.3), and four coefficients whose
        # XSTAR values are less accurate than the others: the Saha factor of type 57 (0.35%) and the
        # photon-density, Milne and Saha-8π coefficients (0.1 to 0.2%)
        looser = Dict(:proton_mass => 3e-2, :eightpi => 1e-2, :saha_coeff_cb => 5e-3, :saha_ci => 6e-3, :fo_saha => 3e-3,
            :bb_coeff => 2e-3, :milne_coeff => 2e-3)
        for name in fieldnames(Radix.Constants)
            @test getfield(U, name) ≈ getfield(K, name) rtol=get(looser, name, 1e-3)
        end
        # fo_saha and bb_coeff only occur together in the recombination term, and their product
        # agrees better than either does
        @test U.fo_saha*U.bb_coeff ≈ K.fo_saha*K.bb_coeff rtol=1.5e-3
        # ... and not identical: it is a different set
        @test U != K
        @test U.kT_eV == Float64(0.861707f0)   # ucalc's single-precision literal
    end

    @testset "hydrogenic decay rates" begin
        # NIST: A(2p→1s) = 6.2649e8, A(3p→2s) = 2.2449e7, A(3p→1s) = 1.6725e8 s⁻¹ for hydrogen; He⁺ scales as Z⁴
        # (the coefficient has the reduced mass of hydrogen, as XSTAR's does)
        @test Radix.anl1(2, 1, 0, 1)[2] ≈ 6.2649e8 rtol=1e-4
        @test Radix.anl1(3, 2, 0, 1)[2] ≈ 2.2449e7 rtol=1e-4
        @test Radix.anl1(3, 1, 0, 1)[2] ≈ 1.6725e8 rtol=1e-4
        @test Radix.anl1(2, 1, 0, 2)[2] ≈ 16*6.2649e8 rtol=1e-4
    end

    @testset "recombination coefficients agree" begin
        # The spontaneous recombination integral of the photoionization code (Saha factor, 8π, photon
        # density 1/(h³c²)) and the Milne integral (4π/((2π m)^{3/2} c²)) are two derivations of the same
        # coefficient from the same cross section; for a Kramers cross section on a database-like grid
        # they agree to the accuracy of the integrations (1-2%), which would expose a wrong factor.
        E = [0.1*exp(0.0015*(i - 1)) for i in 1:9999]
        rad = Radix.Radiation(E, zeros(length(E)))
        Eth, Ry = 13.6, K.Ry_eV
        ε = vcat(0.0, 0.003 .* 1.35 .^ (0:30))                # Ry above the threshold
        σ = 6.3e-18 .* (Eth ./ (Eth .+ ε*Ry)).^3              # cm²
        for T in (0.5, 1.0, 2.0)                              # 10⁴ K
            swrat = 0.5
            r = Radix.photoionization_integrals_hunt(rad, Eth, ε, σ, T, swrat, 1.0)
            α = Radix.milne_recombination(T*1e4, ε, σ/Radix.Mb, Eth/Ry)
            @test r.rrrt/swrat ≈ α rtol=3e-2
        end
    end

    @testset "other CODATA sets" begin
        K2018 = Radix.Constants(CODATA2018)
        @test K2018.kB_cgs == K.kB_cgs                  # exact since the 2019 SI
        @test K2018.Ry_eV ≈ K.Ry_eV rtol=1e-8
        @test K2018.Ry_eV != K.Ry_eV
    end

    @testset "copying with changes" begin
        K3 = Radix.Constants(K; kT_eV=0.5)
        @test K3.kT_eV == 0.5 && K3.Ry_eV == K.Ry_eV
        @test Radix.Constants(NamedTuple{fieldnames(Radix.Constants)}(getfield.(Ref(K), fieldnames(Radix.Constants)))) == K
    end

    @testset "switching the constants" begin
        @test Radix.with_constants(() -> Radix.constants().kT_eV, U) == Float64(0.861707f0)
        @test Radix.constants() === K                     # restored
        @test_throws ErrorException Radix.with_constants(() -> error("x"), U)
        @test Radix.constants() === K                     # restored after an error
    end

    @testset "the rates use them" begin
        # ElectronImpact2 against ucalc's own values: close but not identical with CODATA
        ndiff = 0
        for (c, _) in Iterators.take(ucalc_cases("type98"), 60)
            coef, cell, lv = ucalc_inputs(c)
            r = Radix.rate(coef, cell; levels=lv)
            u = Radix.with_constants(() -> Radix.rate(coef, cell; levels=lv), U)
            r.frate == 0 && continue
            # (the Boltzmann factor amplifies the 3e-5 difference in k by ΔE/kT)
            @test r.frate ≈ u.frate rtol=1e-2
            @test r.irate ≈ u.irate rtol=1e-2
            ndiff += r.frate != u.frate
        end
        @test ndiff > 0
        # the photoionization integrals use 4π and 8π: the photoionization rate cancels them, the
        # recombination does not (ucalc's 8π is 0.7% high)
        nrec = 0
        for (c, _) in Iterators.take(ucalc_cases("type59"), 80)
            coef, cell, lv = ucalc_inputs(c)
            call() = Radix.rate(coef, cell; levels=lv, radiation=ucalc_radiation(c), nlev=c.nlev, abund=(c.cond[9], c.cond[10]))
            r = call()
            u = Radix.with_constants(call, Radix.Constants(Radix.constants(); fourpi=U.fourpi, eightpi=U.eightpi))
            r.frate == 0 && continue
            @test r.frate ≈ u.frate rtol=1e-12
            r.irate == 0 || @test r.irate ≈ u.irate rtol=1e-2
            nrec += r.irate != u.irate
        end
        @test nrec > 0
        # detailed balance holds with either set
        for (c, _) in Iterators.take(ucalc_cases("type56"), 40)
            coef, cell, lv = ucalc_inputs(c)
            r = Radix.rate(coef, cell; levels=lv)
            r.init == 0 && continue
            lo, up = lv[(coef.ion, r.init)], lv[(coef.ion, r.final)]
            expected = Float64(up.g)/Float64(lo.g)*Radix.expo(-abs(Float64(up.E) - Float64(lo.E))/(Radix.constants().kT_eV*cell.T))
            r.irate == 0 || @test r.frate/r.irate ≈ expected rtol=1e-9
        end
    end
end

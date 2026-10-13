# The energy exchanged in Compton scattering: Klein-Nishina scattering by a Maxwell-Jüttner gas (`thermal_exchange`), the surrogate model of it (the default) and XSTAR's
# interpolation in the table `coheat.dat` (`compton_surrogate = false`, what `ucalc()` has).

@testset "Compton scattering" begin
    @testset "Gauss-Legendre rule" begin
        x, w = Radix.gauss_legendre(8)
        @test sum(w) ≈ 2 rtol=1e-14
        @test sum(w .* x.^14) ≈ 2/15 rtol=1e-13               # exact up to degree 15
        @test sum(w .* x.^3) ≈ 0 atol=1e-14
    end

    @testset "thermal exchange from first principles" begin
        # the Thomson limit: ⟨Δε⟩/ε = 4 Te - E (a cold gas gains 4 kT per photon energy and loses the recoil E)
        for (Te, E) in ((1e-6, 1e-6), (1e-5, 1e-7), (1e-4, 1e-6), (1e-7, 1e-5))
            @test Radix.thermal_exchange(Te, E) ≈ 4Te - E rtol=2e-3
        end
        # the recoil of a cold gas is the Klein-Nishina one: the gas takes energy from a hard photon
        @test Radix.thermal_exchange(1e-6, 1.0) < 0
        # an ultrarelativistic gas gives a soft photon 4/3 ⟨γ² β²⟩ = 16 Te² of its energy (⟨γ²⟩ = 12 Te²), not 4 Te
        @test Radix.thermal_exchange(10.0, 1e-6) ≈ 16*10.0^2 rtol=0.2
        # converged: more points change it by less than the tolerance of the surrogate
        for (Te, E) in ((1e-3, 1e-3), (0.02, 3e-4), (1.0, 0.1))
            @test Radix.thermal_exchange(Te, E; n=64) ≈ Radix.thermal_exchange(Te, E; n=96) rtol=3e-4
        end
    end

    @testset "the surrogate follows it" begin
        K = Radix.constants()
        @test K.compton_surrogate && !Radix.ucalc().compton_surrogate
        compton = Radix.Compton([1e-7, 1e2], [1e-7, 1e2], zeros(2, 2))
        scale(Te, E, ex) = max(abs(ex), 4Te + E)
        worst = 0.0
        for Te in 10.0 .^ range(-8.5, 0.4; length=11), E in 10.0 .^ range(-7.5, 0.4; length=11)
            ex = Radix.thermal_exchange(Te, E)
            worst = max(worst, abs(Radix.σ(compton, E, Te) - ex)/scale(Te, E, ex))
        end
        @test worst < 1e-3                              # the table is wrong by 90% of that scale in places
        # below and above the range of the surrogate
        @test Radix.σ(compton, 1e-10, 1e-12) ≈ 4e-12 - 1e-10 rtol=1e-2          # 4 Te - E (the edge of the range gives 3e-5 of its scale)
        @test Radix.σ(compton, 1e-5, 1e-11) ≈ 4e-11 - 1e-5 rtol=1e-3
        @test Radix.σ(compton, 50.0, 1e-3) == Radix.constants(() -> Radix.σ(compton, 50.0, 1e-3), Radix.Constants(K; compton_surrogate=false))     # (the table, which has no values here: zeros)
        # differentiable
        f(Te) = Radix.σ(compton, 0.01, Te)
        h = 1e-6*0.05
        @test ForwardDiff.derivative(f, 0.05) ≈ (f(0.05 + h) - f(0.05 - h))/(2h) rtol=1e-4
    end

    @testset "the coefficients are saved, not trained at each start" begin
        # the file of the repository is what a new training gives (to the rounding of the linear solve)
        @test isfile(Radix.compton_surrogate_file)
        saved = Radix.load_compton_surrogate()
        trained = Radix.train_compton_surrogate()
        for a in range(-18.0, 2.0; length=7), b in range(-16.0, 2.0; length=7)
            @test saved((a, b)) ≈ trained((a, b)) rtol=1e-6 atol=1e-9
        end
        @test saved.coeff ≈ trained.coeff rtol=1e-6
        # written and read back
        mktemp() do path, io
            close(io)
            Radix.save_compton_surrogate(path; nodes=(11, 10))
            small = Radix.load_compton_surrogate(path)
            @test length(small.x) == 110
            @test small((-5.0, -3.0)) ≈ Radix.train_compton_surrogate(; nodes=(11, 10))((-5.0, -3.0)) rtol=1e-9
        end
        # what the model uses is the one of the file
        @test Radix.compton_surrogate()((-5.0, -3.0)) == saved((-5.0, -3.0))
    end

    @testset "the table is the option" begin
        coheat = Radix.load(Radix.Compton, joinpath(@__DIR__, "..", "data", "coheat.dat"))
        @test size(coheat.Es) == (101, 101)
        Ktable = Radix.Constants(Radix.constants(); compton_surrogate=false)
        # the table path reproduces the table at its nodes, and the surrogate is within 1% of it at most of them (it is good to 0.1% in the typical case)
        close = 0
        for (i, j) in ((30, 30), (60, 60), (80, 80), (40, 90), (70, 60))
            Te, E = coheat.Te[i], coheat.Eph[j]
            @test Radix.constants(() -> Radix.σ(coheat, E, Te), Ktable) ≈ coheat.Es[i, j] rtol=1e-2
            close += isapprox(Radix.σ(coheat, E, Te), coheat.Es[i, j]; rtol=1e-2)
        end
        @test close >= 4
        # a relativistic gas and a soft photon: the table has the limit 4 Te - E below 10⁻⁴, the surrogate the exchange
        Te, E = 1.0, 1e-6
        @test Radix.constants(() -> Radix.σ(coheat, E, Te), Ktable) ≈ 4Te - E
        @test Radix.σ(coheat, E, Te) > 3*(4Te - E)
        # the cooling of a power law at 10⁸ K: the table and the surrogate differ by a few percent for a soft spectrum
        E = Radix.xstar_energy_grid(999)
        rad = Radix.point_source(E, Radix.powerlaw(E, -2.0, 1e38), 1e12)
        cool_s = Radix.integral(coheat, rad, 1e4)[2]
        cool_t = Radix.constants(() -> Radix.integral(coheat, rad, 1e4)[2], Ktable)
        @test cool_s ≠ cool_t
        @test cool_s ≈ cool_t rtol=0.1
    end
end

# Tests for the QuantDiffForecast-compatible layer (error_models.jl,
# growth_models.jl, odefit.jl, ensemble.jl).
#
#   julia --project=. -e 'using Pkg; Pkg.test()'
# or
#   julia --project=. test/runtests.jl
#
# Reference objective values were computed independently with
# scipy.special.gammaln using QDF's formulas.

using Test
using Random
using Statistics
using DelimitedFiles

isdefined(Main, :MechFit) || include(joinpath(@__DIR__, "..", "src", "MechFit.jl"))
using .MechFit

@testset "QDF-compatible layer" begin

@testset "error_objective matches QDF formulas" begin
    mu = [2.0, 5.0, 10.0]
    y = [3.0, 4.0, 12.0]
    @test error_objective(:normal, mu, y) ≈ 1 + 1 + 4
    @test error_objective(:laplace, mu, y) ≈ 1 + 1 + 2
    @test error_objective(:poisson, mu, y) ≈ -19.148214307344787 rtol = 1e-10
    @test error_objective(:nb1, mu, y; alpha=0.5) ≈ -18.630720321568138 rtol = 1e-10
    @test error_objective(:nb2, mu, y; alpha=0.5) ≈ -17.429181019994534 rtol = 1e-10
    @test error_objective(:nbd, mu, y; alpha=0.5, d=1.5) ≈ -18.17250158577416 rtol = 1e-10
    # NB with d = 1 is NB1
    @test error_objective(:nbd, mu, y; alpha=0.5, d=1.0) ≈ error_objective(:nb1, mu, y; alpha=0.5)
    # zeros in the curve become 0.001; invalid curves give Inf
    @test error_objective(:poisson, [0.0, 1.0], [0.0, 1.0]) ≈ -(0 * log(0.001) - 0.001) - (0 - 1)
    @test error_objective(:poisson, [0.0, 0.0], [1.0, 1.0]) == Inf
    @test error_objective(:poisson, [NaN, 1.0], [1.0, 1.0]) == Inf
    @test error_objective(:poisson, [-1.0, 1.0], [1.0, 1.0]) == Inf
end

@testset "parameter counts and AICc" begin
    @test n_noise_params(:poisson) == 0
    @test n_noise_params(:nb1) == 1
    @test n_noise_params(:nbd) == 2
    @test n_estimated(:poisson, 1, 0) == 1
    @test n_estimated(:normal, 1, 0) == 2
    @test n_estimated(:nb1, 2, 1) == 4
    @test n_estimated(:nbd, 2, 0) == 4
    n, k = 17, 2
    corr = 2k + 2k * (k + 1) / (n - k - 1)
    a = aicc_qdf(:poisson, -100.0, k, n)
    @test a[1] ≈ -200 + corr
    @test a[1] ≈ a[2] + a[3]
    a = aicc_qdf(:normal, 50.0, k, n)
    @test a[1] ≈ n * log(50 / n) + n * log(2π) + n + corr
    a = aicc_qdf(:laplace, 20.0, k, n)
    @test a[1] ≈ 2n * log(2 * 20 / n) + 2n + corr
    @test_throws ErrorException aicc_qdf(:poisson, 1.0, 5, 6)
    @test qdf_method_code(:nb1) == 3
    @test_throws ErrorException check_error_model(:gaussian)
end

@testset "observation curve and quantiles" begin
    S = reshape([5.0, 7.0, 12.0], :, 1)
    @test vec(observation_curve(S, [1], [true])) == [5.0, 2.0, 5.0]
    @test vec(observation_curve(S, [1], [false])) == [5.0, 7.0, 12.0]
    S2 = hcat([1.0, 2.0, 3.0], [10.0, 15.0, 30.0])
    @test observation_curve(S2, [2, 1], [true, false]) == hcat([10.0, 5.0, 15.0], [1.0, 2.0, 3.0])
    # MATLAB rule: knots at (k - 0.5)/n
    @test matlab_quantile([1.0, 2.0, 3.0, 4.0], 0.3) ≈ 1.7
    @test matlab_quantile([4.0, 1.0, 3.0, 2.0], 0.5) ≈ 2.5
    @test matlab_quantile([1.0, 2.0, 3.0, 4.0], 0.01) == 1.0
    @test matlab_quantile([1.0, 2.0, 3.0, 4.0], 0.99) == 4.0
end

@testset "noise" begin
    rng = Xoshiro(1)
    c = fill(50.0, 4001)
    c[1] = 7.0
    for dist in (:normal, :poisson, :nb1, :nb2, :nbd, :laplace, :nbfactor)
        x = add_noise(c, dist, dist == :normal ? 3.0 : 0.5, 1.5, rng)
        @test x[1] == 7.0
        @test abs(mean(x[2:end]) - 50) < 2
    end
    x = add_noise(c, :nb1, 2.0, 1.0, rng)[2:end]
    @test var(x) ≈ 50 + 2 * 50 rtol = 0.15
    @test variance_mean_ratio(fill(5.0, 30)) == 1.0     # floored at 1
    @test nb_draw(rng, 3.0, 3.0) isa Integer             # falls back to Poisson
end

@testset "WIS and performance" begin
    s = collect(1.0:100.0)
    # an observation at the median scores lower than one far outside
    @test wis_qdf(s, 50.5) < wis_qdf(s, 200.0)
    draws = repeat(collect(1.0:5.0), 1, 200) .+ 0.0
    perf = ode_performance(draws, collect(1.0:5.0), 3, 2)
    @test perf.calibration.MAE == 0
    @test perf.calibration.coverage == 100
    @test length(perf.forecast) == 2
    @test perf.forecast[2].horizon == 2
    @test isempty(ode_performance(draws, collect(1.0:4.0), 3, 2).forecast)   # obs too short
    q = draw_quantiles(draws)
    @test q.median == collect(1.0:5.0)
    @test size(q.quantiles) == (5, length(QDF_QUANTILE_LEVELS))
end

@testset "ensemble weights, pooling, baselines" begin
    w = ensemble_weights(:equal, [10.0, 20.0, 30.0], [1.0, 2.0, 3.0])
    @test w ≈ fill(1 / 3, 3)
    for rule in (:akaike, :inverse_wis, :softmax_wis)
        w = ensemble_weights(rule, [10.0, 12.0, 30.0], [1.0, 2.0, 3.0])
        @test sum(w) ≈ 1
        @test issorted(w; rev=true)          # better score, larger weight
    end
    @test ensemble_weights(:akaike, [10.0, 12.0], [1.0, 1.0]) ≈ [1, exp(-1)] ./ (1 + exp(-1))
    @test_throws ErrorException ensemble_weights(:median, [1.0], [1.0])

    members = [fill(Float64(k), 4, 50) for k in 1:3]
    P = linear_pool([0.5, 0.3, 0.2], members; n_ens=10, seed=1)
    @test size(P) == (4, 10)
    @test count(==(1.0), P[1, :]) == 5
    @test count(==(2.0), P[1, :]) == 3
    @test count(==(3.0), P[1, :]) == 2

    y = [100.0, 120.0, 110.0, 130.0, 125.0, 140.0, 150.0]
    F = baseline_forecast(y, 4, :flat; n_draws=2000, seed=1)
    @test size(F) == (11, 2000)
    @test F[1:7, 1] == y
    @test matlab_quantile(F[9, :], 0.5) ≈ y[end] atol = 1e-8
    T = baseline_forecast(y, 4, :trend; n_draws=2000, seed=1)
    @test matlab_quantile(T[11, :], 0.5) > y[end]          # upward drift carried forward
    @test all(T .>= 0)
end

@testset "growth model right-hand sides" begin
    du = zeros(1)
    ggm!(du, [100.0], [0.5, 0.8], 0.0);  @test du[1] ≈ 0.5 * 100^0.8
    glm!(du, [100.0], [0.5, 1.0, 1000.0], 0.0);  @test du[1] ≈ 0.5 * 100 * 0.9
    rich!(du, [100.0], [0.5, 2.0, 1000.0], 0.0); @test du[1] ≈ 0.5 * 100 * (1 - 0.01)
    grm!(du, [100.0], [0.5, 1.0, 1.0, 1000.0], 0.0); @test du[1] ≈ 0.5 * 100 * 0.9
    lm!(du, [100.0], [0.5, 1000.0], 0.0); @test du[1] ≈ 45.0
    gom!(du, [100.0], [0.5, 0.1], 2.0);   @test du[1] ≈ 50 * exp(-0.2)
    ggm!(du, [-1.0], [0.5, 0.8], 0.0);    @test du[1] == 0       # clipped, no DomainError
    @test r0_seir1([0.6, 0.5, 0.25, 1e5]) ≈ 2.4
end

@testset "synthetic recovery: GGM, Poisson" begin
    truth = FitSpec(model=GGM_MODEL, initial=[0.4, 0.8], lower=[0.0, 0.0], upper=[10.0, 1.0],
                    u0=[5.0], error_model=:poisson)
    t = collect(0.0:24.0)
    z = [0.4, 0.8, 5.0, 0.0, 0.0]
    y = vec(simulate_observation(truth, z, t))
    spec = FitSpec(model=GGM_MODEL, initial=[0.2, 0.5], lower=[0.0, 0.0], upper=[10.0, 1.0],
                   u0=[5.0], error_model=:poisson)
    fit = fit_ode(spec, t, y; n_starts=5, threaded=false)
    @test fit.params ≈ [0.4, 0.8] rtol = 0.01
    @test fit.fitted[1, 1] ≈ y[1]                # fix_I0 anchors the first point
    @test fit.k == 2
    @test isfinite(fit.aicc)
end

@testset "synthetic recovery: SEIR beta, NB1" begin
    N = 1e5
    mk(beta) = FitSpec(model=SEIR1_MODEL, initial=[beta, 0.5, 0.25, N], lower=[0.0, 0.0, 0.0, 20.0],
                       upper=[10.0, 2.0, 2.0, 1e6], fixed=[false, true, true, true],
                       u0=[N - 4, 0.0, 4.0, 0.0, 4.0], error_model=:nb1, composite=r0_seir1)
    t = collect(0.0:29.0)
    z = [0.75, 0.5, 0.25, N, 4.0, 0.0, 0.0]
    y = vec(simulate_observation(mk(0.75), z, t))
    fit = fit_ode(mk(0.4), t, y; n_starts=5, threaded=false)
    @test fit.params[1] ≈ 0.75 rtol = 0.01
    @test fit.composite ≈ 3.0 rtol = 0.01
    @test fit.alpha < 0.1                        # noise-free data: little overdispersion
    @test fit.k == 2                             # beta + alpha
end

@testset "matches MATLAB QuantDiffForecast (flu1918, NB1, days 1-17)" begin
    # Reference values from a QDF run (R2025b, validation/qdf_reference_flu1918_nb1.csv).
    # Point estimates and AICc only: bootstrap quantities depend on the RNG.
    raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-flu1918SF.txt"))
    t, y = raw[1:17, 1], raw[1:17, 2]
    I0 = y[1]
    N = 550_000.0
    seir = FitSpec(model=SEIR1_MODEL, initial=[0.6, 1/1.9, 1/4.1, N], lower=[0.0, 0.0, 0.0, 20.0],
                   upper=[10.0, 2.0, 2.0, 1e6], fixed=[false, true, true, true],
                   u0=[N - I0, 0.0, I0, 0.0, I0], error_model=:nb1, composite=r0_seir1)
    f = fit_ode(seir, t, y; threaded=false)
    @test f.params[1] ≈ 0.771832095238869 rtol = 1e-4
    @test f.composite ≈ 3.16451159047936 rtol = 1e-4
    @test f.aicc ≈ -1688.099 atol = 0.01

    ggm = FitSpec(model=GGM_MODEL, initial=[0.18, 0.9], lower=[0.0, 0.0], upper=[10.0, 1.0],
                  u0=[I0], error_model=:nb1)
    f = fit_ode(ggm, t, y; threaded=false)
    @test f.params ≈ [0.69364207433386, 0.764435000968705] rtol = 2e-3
    @test f.aicc ≈ -1707.266 atol = 0.01

    glm = FitSpec(model=GLM_MODEL, initial=[0.3, 0.9, 3620.0], lower=[0.0, 0.0, 20.0],
                  upper=[10.0, 1.0, 1e6], u0=[I0], error_model=:nb1)
    f = fit_ode(glm, t, y; threaded=false)
    @test f.aicc <= -1703.777 + 0.01          # at least as good as QDF's optimum
    @test f.params[1:2] ≈ [0.693300608901869, 0.76457952497592] rtol = 2e-3
end

@testset "forecast + ensemble smoke test" begin
    spec = FitSpec(model=GGM_MODEL, initial=[0.3, 0.9], lower=[0.0, 0.0], upper=[10.0, 1.0],
                   u0=[5.0], error_model=:poisson)
    rng = Xoshiro(3)
    tt = collect(0.0:24.0)
    mu = vec(simulate_observation(spec, [0.4, 0.8, 5.0, 0.0, 0.0], tt))
    y = [mu[1]; [rand(rng, MechFit.Poisson(m)) for m in mu[2:end]]]
    r1 = forecast_ode(spec, tt, y; window=18, horizon=5, B=10, n_starts=2, progress=false, threaded=false)
    r1b = forecast_ode(spec, tt, y; window=18, horizon=5, B=5, n_starts=2, boot_starts=1, progress=false, threaded=false)
    @test r1b.n_success == 5
    spec2 = FitSpec(model=GLM_MODEL, initial=[0.3, 0.9, 5000.0], lower=[0.0, 0.0, 20.0],
                    upper=[10.0, 1.0, 1e6], u0=[5.0], error_model=:poisson)
    r2 = forecast_ode(spec2, tt, y; window=18, horizon=5, B=10, n_starts=2, progress=false, threaded=false)
    @test size(r1.draws[1]) == (23, 200)
    @test length(r1.performance[1].forecast) == 5
    E = run_ensemble([r1, r2], ["GGM", "GLM"]; rule=:inverse_wis)
    @test sum(E.weights) ≈ 1
    @test size(E.draws, 1) == 23
    @test length(E.relative_wis.ensemble_flat) == 5
    @test all(isfinite, E.relative_wis.ensemble_trend)
end

end

# compare_qdf.jl
#
# Side-by-side check of MechFit against MATLAB QuantDiffForecast (QDF) on the
# 1918 SF flu data: NB1 errors (QDF method1 = 3), calibration days 1-17,
# horizon 10, B = 300. QDF's numbers are in qdf_reference_flu1918_nb1.csv
# (produced by matlab/run_mechfit_comparison.m).
#
#   julia --threads=auto --project=. validation/compare_qdf.jl
#
# Point estimates and AICc should agree closely. Bootstrap-based quantities
# (WIS, coverage) can differ by a few percent because MATLAB and Julia use
# different random number streams.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
isdefined(Main, :MechFit) || include(joinpath(@__DIR__, "..", "src", "MechFit.jl"))
using .MechFit
using DelimitedFiles
using Printf

raw = readdlm(joinpath(@__DIR__, "..", "data", "curve-flu1918SF.txt"))
t_all, y_all = raw[:, 1], raw[:, 2]
I0, N = y_all[1], 550_000.0

specs = [
    "SEIR" => FitSpec(model=SEIR1_MODEL, initial=[0.6, 1/1.9, 1/4.1, N], lower=[0.0, 0.0, 0.0, 20.0],
                      upper=[10.0, 2.0, 2.0, 1e6], fixed=[false, true, true, true],
                      u0=[N - I0, 0.0, I0, 0.0, I0], error_model=:nb1, composite=r0_seir1, composite_name="R0"),
    "GGM"  => FitSpec(model=GGM_MODEL, initial=[0.18, 0.9], lower=[0.0, 0.0], upper=[10.0, 1.0],
                      u0=[I0], error_model=:nb1),
    "GLM"  => FitSpec(model=GLM_MODEL, initial=[0.3, 0.9, 3620.0], lower=[0.0, 0.0, 20.0],
                      upper=[10.0, 1.0, 1e6], u0=[I0], error_model=:nb1),
]

ref = readdlm(joinpath(@__DIR__, "qdf_reference_flu1918_nb1.csv"), ','; header=true)[1]
qdf(model, q) = Float64(ref[findfirst(i -> ref[i, 1] == model && ref[i, 2] == q, axes(ref, 1)), 3])

rows = Any[]
for (name, spec) in specs
    res = forecast_ode(spec, t_all, y_all; window=17, horizon=10, B=300, seed=1, progress=false)
    f = res.fit
    vals = Dict{String,Float64}()
    for (i, lab) in enumerate(spec.model.param_labels)
        spec.fixed[i] || (vals[lab] = f.params[i])
    end
    spec.composite === nothing || (vals["R0"] = f.composite)
    vals["alpha"] = f.alpha
    vals["AICc"] = f.aicc
    c = res.performance[1].calibration
    fc = res.performance[1].forecast[end]
    vals["calibration_WIS"] = c.WIS
    vals["calibration_coverage"] = c.coverage
    vals["forecast_WIS_h10"] = fc.WIS
    vals["forecast_coverage_h10"] = fc.coverage
    for i in axes(ref, 1)
        ref[i, 1] == name || continue
        q = ref[i, 2]
        haskey(vals, q) && push!(rows, (name, q, Float64(ref[i, 3]), vals[q]))
    end
end

println()
@printf("%-5s %-22s %16s %16s %10s\n", "model", "quantity", "QDF (MATLAB)", "MechFit", "rel. diff")
for (m, q, a, b) in rows
    rd = a == 0 ? abs(b) : abs(b - a) / abs(a)
    @printf("%-5s %-22s %16.6g %16.6g %9.2g%%\n", m, q, a, b, 100rd)
end

outdir = joinpath(@__DIR__, "..", "results", "validation")
mkpath(outdir)
open(joinpath(outdir, "qdf_vs_mechfit_flu1918_nb1.csv"), "w") do io
    println(io, "model,quantity,qdf,mechfit")
    for (m, q, a, b) in rows
        println(io, "$m,$q,$a,$b")
    end
end
println("\nWritten to results/validation/qdf_vs_mechfit_flu1918_nb1.csv")

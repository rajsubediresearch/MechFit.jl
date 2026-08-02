# jalisco_data.jl
#
# Loads the Jalisco measles data files (same schema as the Python pipeline's
# inputs_clean/) and builds the `fixed` NamedTuple + observed-cases matrix
# needed by simulate_epidemic_age. Mirrors 04_fit.py / 05_counterfactual.py's
# load_inputs / build functions.

using CSV
using DataFrames
using Dates

const JALISCO_BANDS = ["0-4", "5-9", "10-19", "20-29", "30-49", "50+"]

"""
    load_jalisco_inputs(in_dir; latent_weeks=1.5, infectious_weeks=1.0, vac_eff=0.93,
                         dose_override=nothing)

Loads cases, doses, population, susceptibility, and the contact matrix from
`in_dir`, puts cases and doses on one common weekly clock (union of both,
zero-filled where one series doesn't cover a given week -- e.g. doses
starting before any cases were observed), and builds the `fixed` NamedTuple
for simulate_epidemic_age.

`dose_override`, if given (a n_weeks x n_age matrix on the same weekly
index grid as the loaded data), replaces the loaded dose schedule when
building `fixed` -- used for counterfactual scenarios (delayed campaign,
reduced intensity, leave-one-band-out) without re-reading the CSVs. The
raw loaded `doses_full` is always returned too, for reporting purposes.

Returns a NamedTuple: (N, S0f, S0lo, S0hi, C, obs, doses_full, all_weeks,
fixed, seed_week, n_weeks).
"""
function load_jalisco_inputs(in_dir::AbstractString; latent_weeks::Real=1.5,
                              infectious_weeks::Real=1.0, vac_eff::Real=0.93,
                              dose_override::Union{Nothing,AbstractMatrix}=nothing)
    cases_df = DataFrame(CSV.File(joinpath(in_dir, "01_cases_weekly_by_age.csv")))
    doses_df = DataFrame(CSV.File(joinpath(in_dir, "01_doses_weekly_by_age.csv")))
    pop_df = DataFrame(CSV.File(joinpath(in_dir, "01_population_by_age.csv")))
    sero_df = DataFrame(CSV.File(joinpath(in_dir, "01_baseline_susceptibility_by_age.csv")))
    C_df = DataFrame(CSV.File(joinpath(in_dir, "02_contact_matrix_6band.csv")))

    for b in JALISCO_BANDS
        b in names(cases_df) || error("cases file missing band column '$b'")
        b in names(doses_df) || error("doses file missing band column '$b'")
    end

    # ---- contact matrix: verify row/column order matches JALISCO_BANDS ----
    row_labels = string.(C_df[:, 1])
    row_labels == JALISCO_BANDS ||
        error("contact matrix row order $row_labels != $JALISCO_BANDS -- order " *
              "matters, a mismatch silently transposes mixing")
    col_names = names(C_df)[2:end]
    col_names == JALISCO_BANDS ||
        error("contact matrix column order $col_names != $JALISCO_BANDS")
    C = Matrix{Float64}(C_df[:, 2:end])

    N = [Float64(pop_df[pop_df.age_band .== b, :N][1]) for b in JALISCO_BANDS]
    S0f = [Float64(sero_df[sero_df.age_band .== b, :suscept][1]) for b in JALISCO_BANDS]
    S0lo = [Float64(sero_df[sero_df.age_band .== b, :suscept_lo][1]) for b in JALISCO_BANDS]
    S0hi = [Float64(sero_df[sero_df.age_band .== b, :suscept_hi][1]) for b in JALISCO_BANDS]

    # ---- common weekly clock ----
    to_dates(col) = eltype(col) <: Dates.Date ? collect(col) : Date.(string.(col))
    cases_weeks = to_dates(cases_df.week)
    doses_weeks = to_dates(doses_df.week)
    all_weeks = sort(unique(vcat(cases_weeks, doses_weeks)))

    spacings = unique(Dates.value.(diff(all_weeks)))
    spacings == [7] || error("common clock has row spacing $spacings days, not a " *
                             "uniform 7 -- cases and doses are on different day anchors")

    n_weeks_total = length(all_weeks)
    obs = zeros(Float64, n_weeks_total, length(JALISCO_BANDS))
    doses_full = zeros(Float64, n_weeks_total, length(JALISCO_BANDS))

    cases_idx = Dict(d => i for (i, d) in enumerate(cases_weeks))
    doses_idx = Dict(d => i for (i, d) in enumerate(doses_weeks))
    for (k, wk) in enumerate(all_weeks)
        if haskey(cases_idx, wk)
            obs[k, :] .= [Float64(cases_df[cases_idx[wk], b]) for b in JALISCO_BANDS]
        end
        if haskey(doses_idx, wk)
            doses_full[k, :] .= [Float64(doses_df[doses_idx[wk], b]) for b in JALISCO_BANDS]
        end
    end

    seed_week = Float64(findfirst(==(minimum(cases_weeks)), all_weeks) - 1)  # 0-based index

    vax = VaxSchedule(collect(0.0:(n_weeks_total - 1)),
                       dose_override === nothing ? doses_full : Matrix{Float64}(dose_override))
    fixed = default_fixed_age(N, S0f, C, vax, n_weeks_total - 1;
                               latent_weeks=latent_weeks, infectious_weeks=infectious_weeks,
                               vac_eff=vac_eff)
    fixed = merge(fixed, (seed_week=seed_week,))

    return (N=N, S0f=S0f, S0lo=S0lo, S0hi=S0hi, C=C, obs=obs, doses_full=doses_full,
            all_weeks=all_weeks, fixed=fixed, seed_week=seed_week, n_weeks=n_weeks_total)
end

"""
    build_fixed(D, S0f, dose_matrix)

Build a fresh `fixed` NamedTuple reusing already-loaded data (N, C,
n_weeks, seed_week from a prior load_jalisco_inputs call `D`), varying only
S0f and the dose schedule. Avoids re-reading the CSVs from disk on every
draw in a bootstrap/uncertainty loop -- mirrors 06_uncertainty.py's _fixed().
"""
function build_fixed(D, S0f::AbstractVector, dose_matrix::AbstractMatrix)
    vax = VaxSchedule(collect(0.0:(D.n_weeks - 1)), Matrix{Float64}(dose_matrix))
    fx = default_fixed_age(D.N, S0f, D.C, vax, D.n_weeks - 1)
    return merge(fx, (seed_week=D.seed_week,))
end

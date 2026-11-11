# growth_models.jl
# Phenomenological growth models as one-state ODEs, in the form QuantDiffForecast
# (QDF; Chowell, Bleichrodt & Luo, Stat Med 2024, doi:10.1002/sim.10036) uses
# for its model files. EXP, GGM, GLM and SEIR1 match QDF's EXP.m, GGM.m, GLM.m
# and SEIR1.m at commit fcd7599 (2026-11-09, archived in validation/qdf_snapshot/);
# Richards, GRM, logistic and Gompertz are not in that QDF version and follow the
# standard definitions used in Chowell's growth-model toolboxes. All of them
# run through the same generic fitter (odefit.jl), bootstrap, forecast
# and ensemble code as the mechanistic models. The single state C is
# cumulative cases; fit it with fit_diff = true to match incidence data.
#
# Parameters are passed as a plain vector p in the order listed below
# (QDF's params0 order). Negative C is clipped at 0 before fractional powers
# so a stray solver step cannot raise a DomainError.

"EXP: dC/dt = r C.  p = [r]"
function exp_growth!(du, u, p, t)
    du[1] = p[1] * u[1]
    return nothing
end

"GGM: dC/dt = r C^p.  p = [r, p]"
function ggm!(du, u, p, t)
    C = max(u[1], zero(u[1]))
    du[1] = p[1] * C^p[2]
    return nothing
end

"GLM: dC/dt = r C^p (1 - C/K).  p = [r, p, K]"
function glm!(du, u, p, t)
    C = max(u[1], zero(u[1]))
    du[1] = p[1] * C^p[2] * (1 - C / p[3])
    return nothing
end

"Richards: dC/dt = r C (1 - (C/K)^a).  p = [r, a, K]"
function rich!(du, u, p, t)
    C = max(u[1], zero(u[1]))
    du[1] = p[1] * C * (1 - (C / p[3])^p[2])
    return nothing
end

"Generalized Richards: dC/dt = r C^p (1 - (C/K)^a).  p = [r, p, a, K]"
function grm!(du, u, p, t)
    C = max(u[1], zero(u[1]))
    du[1] = p[1] * C^p[2] * (1 - (C / p[4])^p[3])
    return nothing
end

"Logistic: dC/dt = r C (1 - C/K).  p = [r, K]"
function lm!(du, u, p, t)
    du[1] = p[1] * u[1] * (1 - u[1] / p[2])
    return nothing
end

"Gompertz: dC/dt = r C exp(-a t).  p = [r, a]  (uses absolute time t, as in QDF's GOM.m)"
function gom!(du, u, p, t)
    du[1] = p[1] * u[1] * exp(-p[2] * t)
    return nothing
end

"QDF's SEIR1: states [S, E, I, R, C], p = [beta, kappa, gamma, N]; C counts E→I onsets."
function seir1!(du, u, p, t)
    β, κ, γ, N = p[1], p[2], p[3], p[4]
    S, E, I = u[1], u[2], u[3]
    inf = β * S * I / N
    du[1] = -inf
    du[2] = inf - κ * E
    du[3] = κ * E - γ * I
    du[4] = γ * I
    du[5] = κ * E
    return nothing
end

"R0 = beta/gamma for seir1! parameters (QDF's R0s.m)."
r0_seir1(p) = p[1] / p[3]

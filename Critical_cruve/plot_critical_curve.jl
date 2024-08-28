#=

    In this file we compute and plot the critical curve.

=#

include("../../main.jl")
using HDF5, Plots, .CriticalCurve, BenchmarkTools

# compute critical curve
spin_case = 0.94; i_case = 60; nPointsCritc = Int(1e5);
alpha, beta = CriticalCurve.compute_critical_curve(spin_case, i_case, nPointsCritc);
alphas = vcat(alpha, alpha)
betas = vcat(beta, -beta)

# rotate critical curve
alpha_rot = zeros(2 * nPointsCritc);
beta_rot = zeros(2 * nPointsCritc);
varphi = 90.0; varphi_rad = deg2rad(varphi);
CriticalCurve.rotate_x_y_arrays!(alphas, betas, alpha_rot, beta_rot, varphi_rad)

plot(alphas, betas, label="Critical Curve", color="blue", linewidth=2, size = (500, 500))
plot!(alpha_rot, beta_rot, label="Rotated", color="red", linewidth=2)
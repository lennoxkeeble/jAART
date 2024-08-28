#=

    In this file we compute critical curve diameters and fractional asymmetries of the critical curve in a dense grid of spin and inclination angles.

=#

include("../../main.jl")
using HDF5, Plots, .CriticalCurve
data_path="./Results/Data/";
mkpath(data_path)

num_spins = 20; num_angles = 20;
spins = range(start=0.01, stop=0.99, length=num_spins) |> collect;
angles = range(start=1.0, stop=89.0, length=num_angles) |> collect;

spin_inc_pairs = vec([(a, i) for a in spins, i in angles]);
num_pairs = length(spin_inc_pairs)

nPointsCritc = 1000000;
threaded_spins = zeros(num_pairs);
threaded_angles = zeros(num_pairs);
d_perp = zeros(num_pairs);
d_para = zeros(num_pairs);
asym = zeros(num_pairs);

# extract diameters of critical curve
@time Threads.@threads for j in 1:num_pairs
    spin_case = spin_inc_pairs[j][1];
    i_case = spin_inc_pairs[j][2];
    x_rot = zeros(2 * nPointsCritc);
    y_rot = zeros(2 * nPointsCritc);
    
    alpha, beta = CriticalCurve.compute_critical_curve(spin_case, i_case, nPointsCritc);
    alphas = vcat(alpha, alpha); betas = vcat(beta, -beta);
    
    threaded_spins[j] = spin_case
    threaded_angles[j] = i_case
    d_perp[j] = CriticalCurve.compute_diameter(alphas, betas, x_rot, y_rot, 0.0)
    d_para[j] = CriticalCurve.compute_diameter(alphas, betas, x_rot, y_rot, 90.0)
    asym[j] = 100 * (1 - d_perp[j] / d_para[j])
end

h5open(data_path*"critical_curve_diameters_jl.h5", "w") do file
    file["spins"] = threaded_spins
    file["angles"] = threaded_angles
    file["dperp"] = d_perp
    file["dpara"] = d_para
    file["asym"] = asym
end
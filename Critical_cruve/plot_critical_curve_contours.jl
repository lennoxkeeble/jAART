#=

    In this file we use the (very low resolution) grid computed in the file "compute_critical_curve_grid.jl" to plot contours of the parallel diameter and fractional asymmetry of the critical curve,
    which can used for black hole parameter inference.

=#

include("../../main.jl")
using HDF5, Plots, .CriticalCurve, LaTeXStrings
data_path="./Results/Data/";

# load data
h5open(data_path*"critical_curve_diameters_jl.h5", "r") do file
    global threaded_spins = file["spins"][:]
    global threaded_angles = file["angles"][:]
    global d_perp = file["dperp"][:]
    global d_para = file["dpara"][:]
    global asym = file["asym"][:]
end
num_spins = length(unique(threaded_spins))
num_angles = length(unique(threaded_angles))

# define data structures for the grid
spins = range(start=0.01, stop=0.99, length=num_spins) |> collect;
angles = range(start=1.0, stop=89.0, length=num_angles) |> collect;
dpara_grid = zeros(num_spins, num_angles);
dperp_grid = zeros(num_spins, num_angles);
asym_grid = zeros(num_spins, num_angles);

# manipulate the data into the form required for the plot function (takes a NxM matrix where the columns are the x-values, the rows are the y-values, and the matrix elements are the z-values)
for i in eachindex(spins)
    spin_mask = threaded_spins .== spins[i]
    
    masked_angles = threaded_angles[spin_mask]
    masked_dpara = d_para[spin_mask]
    masked_dperp = d_perp[spin_mask]
    masked_asym = asym[spin_mask]
    
    # now sort by angle
    angle_perm = sortperm(masked_angles)
    sorted_dpara = masked_dpara[angle_perm]
    sorted_dperp = masked_dperp[angle_perm]
    sorted_asym = masked_asym[angle_perm]
    
    # now assign grid values
    dpara_grid[:, i] = sorted_dpara
    dperp_grid[:, i] = sorted_dperp
    asym_grid[:, i] = sorted_asym
end

# plot contours
gr()
spin_inc_contour = plot!(framestyle=:box)
contour(spin_inc_contour, spins, angles, dpara_grid, levels = 5, color=:blue, clabels=true, cbar=false, lw=1)
contour!(spin_inc_contour, spins, angles, asym, levels = 5, color=:red, clabels=true, cbar=false, lw=1)
xlabel!(spin_inc_contour, "Spin, "*L"a")
ylabel!(spin_inc_contour, "Inclination, "*L"i")
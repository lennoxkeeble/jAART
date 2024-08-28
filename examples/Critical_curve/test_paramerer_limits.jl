#=

    In this file we test the limits of the critical curve computation. In particular, how many 0's and 9's we can get for the spin and inclination at their lower and upper bounds. The summary is that we recommend
    using 1e-5 < spin < 0.999999 (six zeroes), and 1 < inclination < 90.0. In this code we used a fixed number of points 1e6 in the analysis, with the reason being that the code is quite fast for this number 
    of points (~30ms), and slows down by a factor of ~10 for 1e7 points. We also note that one can use BigFloat to obtain results for arbitrarily small, but non-zero, spins, and that the code runs for spin = 1.0, but
    a gap appears in the critical curve near the negative x axis.

=#

include("../../main.jl")
using HDF5, Plots, .CriticalCurve, BenchmarkTools, LaTeXStrings, Printf

########## TESTING SPIN ##########
### LOW SPIN ###

# at spin_case = 1e-7 and above, the critical curve looks shape is clearly wrong by eye
spin_case_1 = 1e-6; i_case_1 = 60; nPointsCritc = Int(1e6);
alpha_1, beta_1 = CriticalCurve.compute_critical_curve(spin_case_1, i_case_1, nPointsCritc);

alphas_1 = vcat(alpha_1, alpha_1)
betas_1 = vcat(beta_1, -beta_1)

plot(alphas_1, betas_1, label="Critical Curve", color="blue", linewidth=2, size = (500, 500))

# however, if one uses BigFloat, you can achieve a much smaller level of spin
spin_case_2 = BigFloat(1e-30);
alpha_2, beta_2 = CriticalCurve.compute_critical_curve(spin_case_2, i_case_1, nPointsCritc);

alphas_2 = vcat(alpha_2, alpha_2)
betas_2 = vcat(beta_2, -beta_2)

plot(alphas_1, betas_1, label="a = %$(spin_case_1)", color="blue", linewidth=2, size = (500, 500))
plot!(alphas_2, betas_2, label="a = "*LaTeXString(@sprintf("%.0e", spin_case_2)), color="red", linewidth=2)

### HIGH SPIN ###

# recommend using up to 6 9 nines, e.g., a < 0.999999. The code will run for higher spins, but as a → 1, the left hand side of the critical curve has a single, increasinly large gap (see the plot below)
spin_case_3 = 1-1e-6; i_case_3 = 60; nPointsCritc = Int(1e6);
alpha_3, beta_3 = CriticalCurve.compute_critical_curve(spin_case_3, i_case_3, nPointsCritc, threads=true);
CriticalCurve.nanmin(abs.(beta_1[alpha_1 .< 0]))

alphas_3 = vcat(alpha_3, alpha_3)
betas_3 = vcat(beta_3, -beta_3)

plot(alphas_3, betas_3, label="Critical Curve", color="blue", linewidth=2, size = (500, 500))

# setting a = 1
spin_case_4 = 1.0; i_case_4 = 60; nPointsCritc = Int(1e6);
alpha_4, beta_4 = CriticalCurve.compute_critical_curve(spin_case_4, i_case_4, nPointsCritc);
CriticalCurve.nanmin(beta_4)

alphas_4 = vcat(alpha_4, alpha_4)
betas_4 = vcat(beta_4, -beta_4)

plot(alphas_4, betas_4, label="Critical Curve", color="blue", linewidth=2, size = (500, 500))

# note that using BigFloat does not improve the result, which suggests that unlike at a=0, the computation isn't running into numerical precision errors as a → 1

########## TESTING INCLINATION ##########

### LOW INCLINATION ###

# recommend using i > 1. At inclinations lower inclinations, gaps again start to form in the critical curve near the x axis
spin_case_5 = 0.5; i_case_5 = 1e0; nPointsCritc = Int(1e6);
alpha_5, beta_5 = CriticalCurve.compute_critical_curve(spin_case_5, i_case_5, nPointsCritc);
CriticalCurve.nanmin(beta_5)

alphas_5 = vcat(alpha_5, alpha_5)
betas_5 = vcat(beta_5, -beta_5)

plot(alphas_5, betas_5, label="Critical Curve", color="blue", linewidth=2, size = (500, 500))

### HIGH INCLINATION ###

# critical curve computation works for i = 90.
spin_case_6 = 0.5; i_case_6 = 90; nPointsCritc = Int(1e6);
alpha_6, beta_6 = CriticalCurve.compute_critical_curve(spin_case_6, i_case_6, nPointsCritc);
CriticalCurve.nanmin(beta_6)

alphas_6 = vcat(alpha_6, alpha_6)
betas_6 = vcat(beta_6, -beta_6)

plot(alphas_6, betas_6, label="Critical Curve", color="blue", linewidth=2, size = (500, 500))

# compare low and high inclination for same spin

plot(alphas_5, betas_5, label="i = $(i_case_5)", color="blue", linewidth=2, size = (500, 500))
plot!(alphas_6, betas_6, label="i = $(i_case_6)", color="red")


########## CRITICAL CURVE AT EXTREMA OF SPIN AND INCLINATION ##########

# low spin low inclination
spin_case_7 = 1e-6; i_case_7 = 1.0; nPointsCritc = Int(1e6);
alpha_7, beta_7 = CriticalCurve.compute_critical_curve(spin_case_7, i_case_7, nPointsCritc);
CriticalCurve.nanmin(beta_7)

alphas_7 = vcat(alpha_7, alpha_7)
betas_7 = vcat(beta_7, -beta_7)

plot(alphas_7, betas_7, label="Critical Curve", color="blue", linewidth=2, size = (500, 500))

# low spin high inclination --- using 1e-6 with inclination 90.0 produces a gap in the critical curve by eye -- so recommend using spin > 1e-5 with inclination 90.0
spin_case_8 = 1e-5; i_case_8 = 90.0; nPointsCritc = Int(1e6);
alpha_8, beta_8 = CriticalCurve.compute_critical_curve(spin_case_8, i_case_8, nPointsCritc);
CriticalCurve.nanmin(beta_8)

alphas_8 = vcat(alpha_8, alpha_8)
betas_8 = vcat(beta_8, -beta_8)

plot(alphas_8, betas_8, label="Critical Curve", color="blue", linewidth=2, size = (500, 500))

# high spin low inclination
spin_case_9 = 1-1e-6; i_case_9 = 1.0; nPointsCritc = Int(1e6);
alpha_9, beta_9 = CriticalCurve.compute_critical_curve(spin_case_9, i_case_9, nPointsCritc);
CriticalCurve.nanmin(beta_9)

alphas_9 = vcat(alpha_9, alpha_9)
betas_9 = vcat(beta_9, -beta_9)

plot(alphas_9, betas_9, label="Critical Curve", color="blue", linewidth=2, size = (500, 500))

# high spin high inclination
spin_case_10 = 1-1e-6; i_case_10 = 90.0; nPointsCritc = Int(1e6);
alpha_10, beta_10 = CriticalCurve.compute_critical_curve(spin_case_10, i_case_10, nPointsCritc);
CriticalCurve.nanmin(beta_10)

alphas_10 = vcat(alpha_10, alpha_10)
betas_10 = vcat(beta_10, -beta_10)

plot(alphas_10, betas_10, label="Critical Curve", color="blue", linewidth=2, size = (500, 500))

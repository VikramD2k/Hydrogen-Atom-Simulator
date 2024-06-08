##
include("Functions.jl")
# With the above all the required functions are defined
a = 0.529 * 10^-10;  # Radius of hydrogen atom in meters
rₘᵢₙ = 0.5*a;     # This is the start point of r
rₘₐₓ = 50 * a;
r_steps = 30;
θ_steps = 25;
φ_steps = 25;
Total_steps = r_steps * θ_steps * φ_steps;
n = 4;   # principle quantum number
l = 3;   # azimuthal quantum number
m = 0;   # magnetic quantum number

##

@time prob_density, x, y, z =get_prob_density(n, l, m, rₘᵢₙ, rₘₐₓ, r_steps, θ_steps, φ_steps);
k = 0.2; # k changes density of points in the plot. increasing k increases density
L_,plot_ = plot_probs1(Total_steps, prob_density, x, y, z, n, l, m, k);

##
# Save Plot if you like it as a html file (interractive)

savefig(plot_, "$n$L_.html")
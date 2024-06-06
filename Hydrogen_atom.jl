##
using Symbolics
using Plots
using Statistics
##
a = 0.529*10^-10 # Radius of hydrogen atom in meters
function Legrande(l::Int, cosθ::Num)
    # Compute the Legendre polynomial. Here cosθ is a symbolic variable
    D_cst_l = Differential(cosθ)^l # This is lᵗʰ derivative operator wrt cst
    legrande = (1 / (2^l * factorial(l))) * (D_cst_l(cosθ^2 - 1)^l)
    legrande = expand_derivatives(legrande)
    return legrande
end
#%%
function AssoLegrande(l::Int, m::Int, cst::Float64)
    # Define the symbolic variable
    @variables cosθ
    # Compute the associated Legendre polynomial
    D_cst_m = Differential(cosθ)^m
    asslegrande = (-1)^m * (1 - cosθ^2)^(m/2) * D_cst_m(Legrande(l, cosθ))
    asslegrande = expand_derivatives(asslegrande)
    # Substitute the numeric value of cst into the symbolic expression
    asslegrande = substitute(asslegrande, cosθ => cst)
    return asslegrande
end
#%%
function Laguerre(q::Int,x::Num)
    # Laguerre function
    # Here x is a symbolic variable
    if q==0
        leg=1;
    else
        Dx_q = Differential(x)^q
        leg = exp(x)*Dx_q(exp(-x)*x^q)/factorial(q)
        leg = expand_derivatives(leg)
    end
    return leg
    # this returns Laguerre polynomial
end
#%%
function AssoLaguerre(p ::Int, q ::Int, k::Float64)
    # Associated Laguerre function
    @variables x
    if p == 0
        Dx = Differential(x)
        asso_lag = (-1)^p * (p+q,x)
    else
        Dx_p = Differential(x)^p
        asso_lag = (-1)^p * Dx_p(Laguerre(p+q, x))
    end
    # this returns Associated Laguerre polynomial
    asso_lag = expand_derivatives(asso_lag)
    asso_lag = substitute(asso_lag, x => k)
    return asso_lag
end
#%%
function Radial(n::Int, l::Int, r::Float64)
    # Define the symbolic variables
    @variables ρ
    # Compute the radial part of the wavefunction
    radial = sqrt((2/(n*a))^3 * factorial(n-l-1) / (factorial(n+l) * 2 * n)) * exp(-r/(n*a)) * (2*r/(n*a))^l
    radial = radial * AssoLaguerre(2*l+1, n-l-1, r)
    radial = expand_derivatives(radial)
    # Substitute temp = 2*r/(n*a)
    ρ_val = 2*r/(n*a)
    
    # Substitute the numeric value of r into the symbolic expression
    radial = substitute(radial,Dict([ρ => ρ_val]))
    
    return radial
end
#%%
function Angular(l::Int, m::Int, cst::Float64, phi::Float64)
    i = Complex{Float64}(0, 1)  # Define the imaginary unit
    angular = sqrt((2*l+1) * factorial(l-m) / (4*pi * factorial(l+m))) * AssoLegrande(l, m, cst) * exp(i * m * phi)
    return angular
end
##
@variables ρ, cosθ, φ # ρ is radial distance, φ is azimuthal angle, θ is the 
# n, l, m are the principle, azimuthal, and magnetic quantum numbers
##% Initialisations
rₘᵢₙ = 0.5* a
rₘₐₓ = 50* a
r_steps = 10
θ_steps = 10
φ_steps = 10
Total_steps = r_steps * θ_steps * φ_steps
n = 1
l = 0
m = 0
x = zeros(Total_steps)
y = x
z = x
#ψ(n,l,m,θ,r,ϕ) = Radial(n, l, r)*Angular(l, m, cos(θ), ϕ)
# ψ = Radial(...) * Angular(...)
# The Radial(...) calls the Associated Laguerre function which in turn calls the Laguerre function
# Simillarly the Radial(...) calls the Associated Legendre function which in turn calls the Legendre function
ψ² = zeros(Total_steps) # Initializing place holder for |ψ|²
global  counter = 1
##
for r = range(rₘᵢₙ, rₘₐₓ, 100)
    for θ = range(0, π, 30)
        for ϕ=range(0, 2*π, 30)
            print(ϕ)
            # substitute and save the values of psi for all
            # different psi, theta and r
            ψ = Radial(n, l, r)*Angular(l, m, cos(θ), ϕ)
            ψ²[counter] = abs(ψ)^2
            # this is the probability of finding an electron of the orbital (n,l,m) at the position (r,θ,ϕ)
            # transforming spherical coordinates to cartesian coordinates and storing them
            x[counter] = r*sin(θ)*cos(ϕ)
            y[counter] = r*sin(θ)*sin(ϕ)
            z[counter] = r*cos(θ)
            counter = counter +1
        end
    end
end
##
ψ²_sum =sum(ψ²) # This should be close to 1
ψ² = ψ²/ψ²_sum # Now it is normalized
# Initialize the points that need to be plotted
X = []
Y = []
Z = []
c=0
for b=1:PH
    if ψ² < 1.1 && ψ² > 0.9
        c = c + 1
        push!(X,x[c])
        push!(Y,y[c])
        push!(Z,z[c])
            # I am storing all the points at which 0.9<Psi/Mpsi<1.1 
            # Sometimes I changed the range to get better plot
    b=b+1

# Now plotting the stored points
Plots.scatter3d(X, Y, Z, markersize = 2, color = :red)
title!("n= l=2 m=2")
xlabel!("x")
ylabel!("y")
zlabel!("z")
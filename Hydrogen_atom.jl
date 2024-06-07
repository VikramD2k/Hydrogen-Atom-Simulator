##
using Symbolics
using Plots
using Statistics
##
function Legrande(l::Int, cosθ::Num) 
    if l == 0
        legrande = 1
    else
        D_cst_l = Differential(cosθ)^l
        legrande = (1 / (2^l * factorial(l))) * D_cst_l((cosθ^2 - 1)^l)
    end
    return expand_derivatives(legrande)
end
##
function AssoLegrande(l::Int, m::Int, cst::Float64)
    @variables cosθ
    if m == 0
        asslegrande = (-1)^m * (1 - cosθ^2)^(m / 2) * Legrande(l, cosθ)
    else
        D_cst_m = Differential(cosθ)^m
        asslegrande = (-1)^m * (1 - cosθ^2)^(m / 2) * D_cst_m(Legrande(l, cosθ))
    end
    asslegrande = expand_derivatives(asslegrande)
    asslegrande = substitute(asslegrande, cosθ => cst)
    return asslegrande
end
##
function Laguerre(q::Int, x::Num)
    if q == 0
        leg = 1
    else
        Dx_q = Differential(x)^q
        leg = exp(x) * Dx_q(exp(-x) * x^q) / factorial(q)
    end
    return expand_derivatives(leg)
end
##
function AssoLaguerre(p::Int, q::Int, k::Float64)
    @variables ρ
    if p == 0
        asso_lag = (-1)^p * Laguerre(p + q, ρ)
    else
        Dρ_p = Differential(ρ)^p
        asso_lag = (-1)^p * Dρ_p(Laguerre(p + q, ρ))
    end
    asso_lag = expand_derivatives(asso_lag)
    asso_lag = substitute(asso_lag, ρ => k)
    return asso_lag
end
##
function Radial(n::Int, l::Int, r::Float64)
    a = 0.529 * 10^-10  # Radius of hydrogen atom in meters
    r_val = 2 * r / (n * a)
    radial = sqrt((2 / (n * a))^3 * factorial(n - l - 1) / (factorial(n + l) * 2 * n)) * exp(-r / (n * a)) * (2 * r / (n * a))^l
    radial = radial * AssoLaguerre(2 * l + 1, n - l - 1, r_val)
    radial = expand_derivatives(radial)
    return radial
end
##
function Angular(l::Int, m::Int, cst::Float64, phi::Float64)
    i = Complex{Float64}(0, 1)
    angular = sqrt((2 * l + 1) * factorial(l - m) / (4 * π * factorial(l + m))) * AssoLegrande(l, m, cst) * exp(i * m * phi)
    return angular
end
##
function get_prob_density(n, l, m, rₘᵢₙ, rₘₐₓ, r_steps, θ_steps, φ_steps)
    x = zeros(Total_steps)
    y = zeros(Total_steps)
    z = zeros(Total_steps)
    ψ² = zeros(Total_steps)  # Initializing placeholder for |ψ|²
    counter = 1
    for r in range(rₘᵢₙ, rₘₐₓ, length=r_steps)
        for θ in range(0, π, length=θ_steps)
            for ϕ in range(0, 2 * π, length=φ_steps)
                ψ = eval(Radial(n, l, r) * Angular(l, m, cos(θ), ϕ)) # Here we got probability amplitudes at the point (r,θ,φ)
                print(ψ, typeof(ψ))
                ψ²[counter] = abs2(ψ)
                x[counter] = r * sin(θ) * cos(ϕ)
                y[counter] = r * sin(θ) * sin(ϕ)
                z[counter] = r * cos(θ)
                counter = counter + 1
            end
        end
    end
    ψ²_sum = sum(ψ²)  # This should be close to 1
    ψ² = ψ² / ψ²_sum  # Now it is normalized
    return ψ², x, y, z
end
##
function plot_probs(Total_steps::Int, ψ², x, y, z)
    X, Y, Z = [], [], []
    c = 0

    for i in 1:Total_steps
        if ψ²[i] < 1.1 && ψ²[i] > 0.9
            c += 1
            push!(X, x[i])
            push!(Y, y[i])
            push!(Z, z[i])
        end
    end
    Plots.scatter3d(X, Y, Z, markersize=2, color=:red)
    title = "n is $n, l is $l, m is $m"
    title!(title)
    xlabel!("x")
    ylabel!("y")
    zlabel!("z")
end
##
# With the above all the required functions are defined
a = 0.529 * 10^-10  # Radius of hydrogen atom in meters
rₘᵢₙ = 0.5 * a
rₘₐₓ = 50 * a
r_steps = 10
θ_steps = 10
φ_steps = 10
Total_steps = r_steps * θ_steps * φ_steps
n = 1   # principle quantum number
l = 0   # azimuthal quantum number
m = 0   # magnetic quantum number
##
prob_density, x, y, z =get_prob_density(n, l, m, rₘᵢₙ, rₘₐₓ, r_steps, θ_steps, φ_steps)
plot_probs(Total_steps, prob_density, x, y, z)
##
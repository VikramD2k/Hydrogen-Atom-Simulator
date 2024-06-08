##
using Symbolics
using Statistics
using BenchmarkTools
using PlotlyJS
using Dierckx
using SparseArrays
##
function Legrande(l::Int, cosθ::Num) 
    """
    Legrande(l::Int, cosθ::Num)

    Compute the Legendre polynomial for a given order `l` and symbolic variable `cosθ`.
    """
    if l == 0
        legrande = 1;
    else
        D_cst_l = Differential(cosθ)^l;
        legrande = (1 / (2^l * factorial(l))) * D_cst_l((cosθ^2 - 1)^l);
    end
    return expand_derivatives(legrande);
end
##
function AssoLegrande(l::Int, m::Int, cst::Float64)
    """
    AssoLegrande(l::Int, m::Int, cst::Float64)

    Compute the Associated Legendre polynomial for a given order `l`, degree `m`, and numeric value at `cst`.
    """
    @variables cosθ
    if m == 0
        asslegrande = Legrande(l, cosθ);
    elseif m > 0
        D_cst_m = Differential(cosθ)^m;
        asslegrande = (-1)^m * (1 - cosθ^2)^(m / 2) * D_cst_m(Legrande(l, cosθ));
    elseif m < 0
        error("The code for orbitals with -ve m values is under progress, please check out other orbitals")
    end
    asslegrande = expand_derivatives(asslegrande);
    asslegrande = substitute(asslegrande, cosθ => cst);
    return asslegrande;
end
##
function Laguerre(q::Int, x::Num)
    """
    Laguerre(q::Int, x::Num)

    Compute the Laguerre polynomial for a given order `q` and symbolic variable `x`.
    """
    if q == 0
        leg = 1;
    else
        Dx_q = Differential(x)^q;
        leg = exp(x) * Dx_q(exp(-x) * x^q) / factorial(q);
    end
    return expand_derivatives(leg);
end
##
function AssoLaguerre(p::Int, q::Int, k::Float64)
    """
    AssoLaguerre(p::Int, q::Int, k::Float64)

    Compute the Associated Laguerre polynomial for a given order `p`, degree `q`, and numeric value at `k`.
    """
    @variables ρ
    if p == 0
        asso_lag = (-1)^p * Laguerre(p + q, ρ);
    else
        Dρ_p = Differential(ρ)^p;
        asso_lag = (-1)^p * Dρ_p(Laguerre(p + q, ρ));
    end
    asso_lag = expand_derivatives(asso_lag);
    asso_lag = substitute(asso_lag, ρ => k);
    return asso_lag;
end
##
function Radial(n::Int, l::Int, r::Float64)
    """
    Radial(n::Int, l::Int, r::Float64)

    Compute the radial part of the wavefunction for given quantum numbers `n`, `l`, and radial distance `r`.
    """
    a = 0.529 * 10^-10;  # Radius of hydrogen atom in meters
    r_val = 2 * r / (n * a);
    radial = sqrt((2 / (n * a))^3 * factorial(n - l - 1) / (factorial(n + l) * 2 * n)) * exp(-r / (n * a)) * (2 * r / (n * a))^l;
    radial = radial * AssoLaguerre(2 * l + 1, n - l - 1, r_val);
    radial = expand_derivatives(radial);
    return radial;
end
##
function Angular(l::Int, m::Int, cst::Float64, phi::Float64)
    """
    Angular(l::Int, m::Int, cst::Float64, phi::Float64)

    Compute the angular part of the wavefunction for given quantum numbers `l`, `m`, numeric value `cst`, and azimuthal angle `phi`.
    """
    i = Complex{Float64}(0, 1);
    angular = sqrt((2 * l + 1) * factorial(l - m) / (4 * π * factorial(l + m))) * AssoLegrande(l, m, cst) * exp(i * m * phi);
    return angular;
end
##
function get_prob_density(n, l, m, rₘᵢₙ, rₘₐₓ, r_steps, θ_steps, φ_steps)
    """
    get_prob_density(n, l, m, rₘᵢₙ, rₘₐₓ, r_steps, θ_steps, φ_steps)

    Compute the probability density `|ψ|²` at Cartesian coordinates for given quantum numbers `n`, `l`, `m`,
        and range parameters for radial distance `r`, polar angle `θ`, and azimuthal angle `φ`.
    """
    if n<=l
        error("n can't be greater than or equal to l")
    end
    if m > l || m < -l
        error("m should lie between -l to l")
    end
    Total_steps = r_steps * θ_steps * φ_steps
    x = SparseArrays.zeros(Total_steps);
    y = SparseArrays.zeros(Total_steps);
    z = SparseArrays.zeros(Total_steps);
    ψ² = SparseArrays.zeros(Total_steps); # Initializing placeholder for |ψ|²
    counter = 1;
    for r in range(rₘᵢₙ, rₘₐₓ, length=r_steps)
        for θ in range(0, π, length=θ_steps)
            for ϕ in range(0, 2 * π, length=φ_steps)
                ψ = (eval(Radial(n, l, r) * Angular(l, m, cos(θ), ϕ))); # Here we got probability amplitudes at the point (r,θ,φ)
                ψ = Symbolics.value(ψ); # Converts {Num} type to {Float64}
                #print("ψ is ", ψ, " Type of ψ is ", typeof(ψ), "\n")
                ψ²[counter] = abs2(ψ);
                x[counter] = r * sin(θ) * cos(ϕ);
                y[counter] = r * sin(θ) * sin(ϕ);
                z[counter] = r * cos(θ);
                counter = counter + 1;
            end
        end
    end
    ψ²_sum = sum(ψ²);  # This should be close to 1
    ψ² = ψ² / ψ²_sum;  # Now it is normalized
    return ψ², x, y, z;
end
##

function plot_probs1(Total_steps::Int, ψ², x, y, z, n, l, m, k)
    X, Y, Z = Float64[], Float64[], Float64[];
    c = 0;
    ψ²_mean = Statistics.mean(ψ²);
    σψ² = Statistics.std(ψ²);
    # k decides how much volume we want to cover in densiy space
    σ⁺ = ψ²_mean + k*σψ²;
    σ⁻ = ψ²_mean - k*σψ²;
    for i in 1:Total_steps
        if ψ²[i] > σ⁺
            c += 1;
            push!(X, x[i]);
            push!(Y, y[i]);
            push!(Z, z[i]);
        end
    end

    scatter3d = PlotlyJS.scatter3d(
        x=X,
        y=Y,
        z=Z,
        mode="markers",
        marker=attr(
            size=3,
            color=ψ²[1:c],
            colorscale="Viridis",
            showscale=true,
            colorbar=attr(
                title="Probability"
                )
            )
    );
    if l==0
        L = "s";
    elseif l == 1
        if m == 0
            L = "pᶻ";
        elseif m == 1
            L = "pˣ";
        else L = "pʸ";
        end
    elseif l == 2
        if m == 0
            L = "dᶻ²";
        elseif m == 1
            L = "dˣᶻ";
        elseif m == -1
            L = "dʸᶻ";
        elseif m == 2
            L = "dˣʸ";
        else
            L = "dˣ²⁻ʸ²";
        end
    elseif l == 3
        L = "f & m = $m";
    elseif l == 4
        L = "g & m = $m";
    elseif l == 5
        L = "h & m =$m";
    else 
        L = "l = $l, m = $m";
    end


    layout = Layout(
        title="Probability Density Scatter Plot of $n$L orbital",
        scene=attr(
            xaxis_title="x",
            yaxis_title="y",
            zaxis_title="z"
        )
    );
    print("Probability Density Scatter Plot of $n$L orbital");
    plot = Plot(scatter3d, layout);
    display(plot);
    return L, plot
end;
    
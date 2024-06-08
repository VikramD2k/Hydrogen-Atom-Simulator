# Hydrogen-Atom-Simulator
This computes and plots the probability density of any orbital of the hydrogen atom \
Feel free to play around with the `sample plots` of the orbitals from `1s` to `4f` \
Just download the `.html` file and open it with any browser \
There is no limit on how high one can go, for example, it is totally possible to get the plot of `n = 10`, `l = 7`, `m = 6`. \
The `Funtions.jl` contains the following functions  
  1. Laguerre: generates Laguerre polynomial of required degree and order 
  2. AssoLaguerre: computes the Associated Laguerre polynomial from Laguerre's output. It also computes the value of this Associated Laguerre polynomial for the given input 
  3. Lagrande: generates Lagrande polynomial of required degree and order 
  4. AssoLagrande: computes the Associated Lagrande polynomial from Lagrande's output. It also computes the value of this Associated Lagrande polynomial for the given input 
  5. Radial: This computes the radial part of wave function calls AssoLagrnade 
  6. Angular: This computes the angular part of wave function calls AssoLaguerre 
  7. get_prob_density: This computes the probability density `|ψ|²` at Cartesian coordinates for given quantum numbers `n`, `l`, `m`, for the radial distance `r`, polar angle `θ`, and azimuthal angle `φ` 
  8. plot_probs1: This plots the the probability densitys and gives a interractive `.html` file that has the 3D scatter plot 
## Limitations of current version of `Funcitons.jl`
In this version, the code can not handle negative values for m, it can handle only non-negative values of m \
Although the plots for `m = +a` and `m = -a` for some `a` are analogus to one another, I feel it is important to mention this limitation 
### If you like my work Do give me a Star 

# Kinetic Assembly Rates Optimizer

This is a tool that optimizes rate constants for an nmer protein of size n.

Dependencies: Catalyst, OrdinaryDiffEq, DiffEqCallbacks, Plots, Optimization, OptimizationOptimisers, NNlib, ReverseDiff
Julia Version: 1.10.0

## Simplified script
To use this code simply run the optimize.jl script with the following arguments:
  
  &ensp;&ensp;Arg1: nmer size (n should be an integer)
  
  &ensp;&ensp;Arg2: maximum number of iterations for optimization (should be an integer)
  
  &ensp;&ensp;Arg3: learning rate (should be a float)
  
  &ensp;&ensp;Arg4: Which automatic differentiator to use ("Reverse","Rev","R","rev","r", or "reverse" for AutoReverseDiff() and anything other string for AutoForwardDiff()
  
  &ensp;&ensp;Arg5: Which integrator to use (Please specify the integrator in julia format. For example if you want to use QNDF, please input "QNDF()")
### Example
Here is an example of how to run the code with a trimer, 1 iteration, .01 learning rate, AutoForwardDiff, and QNDF:
```
julia +1.10.0 ./optimize.jl 3 1 .01 f "QNDF()"
```
This code above automatically sets the initial monomer conentrations to 100 and all forward rates to 10.

## Using your own Julia script for maximum customizability
`ReactionNetwork.jl` contains `get_fc_rn(n::Int)` which will automatically create a rate growth model network of size n based on input.
You can use the code like so:
```
trimer = get_fc_rn(3)
```
`optim.jl` contains `optim(rn::ReactionSystem,tspan::Tuple{Float64,Float64},p_init::Vector{Float64},monomer_conc::Vector{Float64},lr::Float64,iters::Int,AD::AbstractADType,integrator::Any;verbose::Bool=false)`.
- rn: Your reaction network generated above
- tspan: The time span to integrate over to determine yield for calculating loss each iteration. I recommend going until the start of the kinetic trapping or just before kinetic trapping occurs. The longer the time span the longer it will take.
- p_init: Your initial forward rates. Backward rates are calculated based on forward rates. The vector should be of length n-1 where the value at each index i corresponds to the forward rates of species of size i.
- monomer_conc: The initial starting concentrations of each monomer species. The vector should be of length n where the value at each index i corresponds to species Xi.
- lr: The learning rate for Adam optimizer.
- iters: The maximum number of iteration to use for optimization.
- AD: Method of automatic differentiation. Currently only works for `AutoForwardDiff()` and `AutoReverseDiff()`.
- integrator: Input the integrator you wish to use.
- verbose: Whether to print the initial forward and backwards rates.
### Example
Here is an example of running a trimer.
```
new_params = optim(trimer,(0.,.1),[10.,10.],[100.,100.,100.],.01,1000,Optimization.AutoForwardDiff(),KenCarp4())
```
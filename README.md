This is a tool that optimizes rate constants for an nmer protein of size n.

Dependencies: Catalyst, OrdinaryDiffEq, DiffEqCallbacks, Plots, Optimization, OptimizationOptimisers, NNlib, ReverseDiff
Julia Version: 1.10.0

To use this code simply run the trimer.jl script (name to be changed) with the following arguments:
  Arg1 = nmer size (n should be an integer)
  Arg2 = maximum number of iterations for optimization (should be an integer)
  Arg3 = learning rate (should be a float)
  Arg4 = Which automatic differentiator to use ("Reverse","Rev","R","rev","r", or "reverse" for AutoReverseDiff() and anything other string for AutoForwardDiff()
  Arg5 = Which integrator to use (Please specify the integrator in julia format. For example if you want to use QNDF, please input "QNDF()")

Here is an example of how to run the code with a trimer, 1 iteration, .01 learning rate, AutoForwardDiff, and QNDF:
```
julia +1.10.0 ./trimer.jl 3 1 .01 f "QNDF()"
```

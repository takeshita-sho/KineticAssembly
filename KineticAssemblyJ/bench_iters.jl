using Catalyst
using OrdinaryDiffEq, DiffEqCallbacks
using Plots
using BenchmarkTools
using Optimization
include("./optim.jl")
include("./ReactionNetwork.jl")
# Can make this automatically create based on input n
#This defines the system of the homo rate trimer

tspan = (0., .1)
n = parse(Int,ARGS[1])
iters = parse(Int,ARGS[2])
lr= parse(Float64,ARGS[3])
AD = Optimization.AutoForwardDiff()
if ARGS[4] in ["Reverse","Rev","R","rev","r","reverse"]
    AD = Optimization.AutoReverseDiff()
end
integrator = eval(Meta.parse(ARGS[5]))
println(AD)
println("Iters: ",iters)



#for integrator in [Rosenbrock23(),TRBDF2(),QNDF(),Rodas4P(),Rodas5P(),Kvaerno5(),KenCarp4()]
println(integrator)

println("$(n)mer")
println("Peak Memory")
nmer = get_fc_rn(n)
params = fill(10.0,n-1)
monomer_conc = fill(100.0,n)
#optim(nmer,tspan,params,monomer_conc,lr,iters,AD,integrator)
benchmark_result = @btime optim($nmer, $tspan, $params, $monomer_conc, $lr, $iters, $AD,$integrator)
println(benchmark_result)

flush(stdout)
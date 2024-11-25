using Catalyst
using OrdinaryDiffEq, DiffEqCallbacks
using Plots
using BenchmarkTools
using Optimization
include("./optim.jl")
include("./ReactionNetwork.jl")
# Can make this automatically create based on input n
#This defines the system of the homo rate trimer

println(Sys.total_memory())
println(Sys.free_memory())
tspan = (0., .1)
lr=.01
#iters = 10000
AD = Optimization.AutoForwardDiff()
#integrator = QNDF()
#n = 3
println(AD)
iters = 1
for integrator in [Rosenbrock23(),TRBDF2(),QNDF(),FBDF(),Rodas4P(),Rodas5P(),Kvaerno5(),KenCarp4()]
    println(integrator)
    for n in 3:8
        #println("Iterations $(iters)")
        println("$(n)mer")
        nmer = get_fc_rn(n)
        params = fill(10.0,n-1)
        monomer_conc = fill(100.0,n)
        #optim(nmer,tspan,params,monomer_conc,lr,iters,AD,integrator)
        benchmark_result = @btime optim($nmer, $tspan, $params, $monomer_conc, $lr, $iters, $AD,$integrator)
        println(benchmark_result)
        
        flush(stdout)
    end
end
#=
trimer = get_fc_rn(4)


# Define simulation settings
monomer_conc = [100.0,100.0,100.0,100.0] #concentration
tspan = (0., .1) #time span
lr=.01

#These params are for homorates
#deltaG -20
#params = [:k1 => 50.0, :k2 => .0002, :k3 => 50.0, :k4 => 4.24e-12] # initial rates
params = [10.0,10.0,10.0]
iters=10000
println(optim(trimer,tspan,params,monomer_conc,lr,iters))
=#



#=
#Run simulation
ode = ODEProblem(trimer, u0, tspan, params)#; jac = true) #Using the jacobian
# David seemed to use specific sampling to inidicate where to step to
# I might just leave it up to solve to determine
sol = solve(ode,Rodas5P()) # using Rodas5P since better for stiff problems, could also use TRBDF2
println(sol.u[end][end])
println(Array(sol)[end][end])
plot(sol; lw = 5)
=#

#I think Spencer doesnt calculate kinetic traps or time, just plugged in manually, but need to check
#Also I think Davids loss is much more complicated than what Spencer implemented - loss is just negative yield+penalty

#can solve yield at eq if not sure of final yield
#solve system of equations of reactions to get final conce of each species

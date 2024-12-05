using Catalyst
using OrdinaryDiffEq, DiffEqCallbacks
using Plots
using Optimization

using Roots
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
println(integrator)
#for n in 3:10
println("$(n)mer")
nmer = get_fc_rn(n)
println("Generated Reaction Network")
flush(stdout)

params = fill(10.0,n-1)
monomer_conc = fill(100.0,n)
u0 = get_species_conc(monomer_conc,nmer)
#=
ode_init = ODEProblem(nmer, u0, (.00000001,10000), get_rates(params,unique!(reactionrates(nmer))))
sol_init = solve(ode_init,integrator)

plot(sol_init,xaxis=:log; lw = 5,legend=:outerright,title="$(n)mer Init")
savefig("$(n)mer_init.png")
=#

println("Beginning Optimization")
flush(stdout)
new_params = optim(nmer,tspan,params,monomer_conc,lr,iters,AD,integrator)
println("Finished Optimization")
flush(stdout)
#println("New parameters: ",new_params)
#println(new_params)
#flush(stdout)

#maybe just go through all time points


ode = ODEProblem(nmer, u0, (.00000001,10000), new_params)#; jac = true) #Using the jacobian
sol = solve(ode,integrator)
plot(sol,xaxis=:log; lw = 5,legend=:outerright,title="$(n)mer $(iters) Iters")
savefig("$(n)mer$(iters)iters.png")



# David seemed to use specific sampling to inidicate where to step to


#I think Spencer doesnt calculate kinetic traps or time, just plugged in manually, but need to check
#Also Davids loss is mean squared error - Spencers loss is just negative yield+penalty

#can solve yield at eq if not sure of final yield
#solve system of equations of reactions to get final conce of each species

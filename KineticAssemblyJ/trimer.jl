using Catalyst
using DifferentialEquations
using Plots
include("./optim.jl")
include("./ReactionNetwork.jl")
# Can make this automatically create based on input n
#This defines the system of the homo rate trimer
trimer = get_fc_rn(3)


# Define simulation settings
monomer_conc = [100.0,100.0,100.0] #concentration
tspan = (0., .1) #time span
lr=.01

#These params are for homorates
#deltaG -20
#params = [:k1 => 50.0, :k2 => .0002, :k3 => 50.0, :k4 => 4.24e-12] # initial rates
params = [10.0,10.0]
iters=100
println(optim(trimer,tspan,params,monomer_conc,lr,iters))




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

using Catalyst
using DifferentialEquations
using Plots
include("./optim.jl")

# Can make this automatically create based on input n
#This defines the system of the homo rate trimer
trimer = @reaction_network begin
    (k1,k2), A+B <--> AB
    (k1,k2), A+C <--> AC
    (k1,k2), B+C <--> BC
    (k3,k4), AB+C <--> ABC
    (k3,k4), AC+B <--> ABC
    (k3,k4), BC+A <--> ABC   
end


# Define simulation settings
u0 = [:A => 100.0, :B => 100.0, :C => 100.0, :AB => 0.0,
 :BC => 0.0, :AC => 0.0, :ABC => 0.0] #concentration
tspan = (0., .1) #time span

#These params are for homorates
params = [:k1 => 50.0, :k2 => 10.0, :k3 => 50.0, :k4 => 10.0] # initial rates

#print(optim(trimer, tspan,params,u0))





#Run simulation
ode = ODEProblem(trimer, u0, tspan, params)#; jac = true) #Using the jacobian
# David seemed to use specific sampling to inidicate where to step to
# I might just leave it up to solve to determine
sol = solve(ode,Rodas5P()) # using Rodas5P since better for stiff problems, could also use TRBDF2
println(sol.u[end][end])
println(Array(sol)[end][end])
plot(sol; lw = 5)


#I think Spencer doesnt calculate kinetic traps or time, just plugged in manually, but need to check
#Also I think Davids loss is much more complicated than what Spencer implemented - loss is just negative yield+penalty

#can solve yield at eq if not sure of final yield
#solve system of equations of reactions to get final conce of each species

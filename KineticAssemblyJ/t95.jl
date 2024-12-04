using Catalyst
using OrdinaryDiffEq, DiffEqCallbacks
using Plots
using Optimization
using Roots
include("./optim.jl")
include("./ReactionNetwork.jl")
function t95(sol,params)
    for i in 1:length(sol.t)
        if sol.u[i][end] > 95
            m = (sol.u[i][end]-sol.u[i-1][end])/(sol.t[i]-sol.t[i-1])
            b = sol.u[i][end]-m*sol.t[i]
            t = ((95-b)/m)
            kfmax = maximum(values(params))
            #ask adip about this scale time?
            #100 is Cmono - concentration I would like to reach
            return t*kfmax*100
        end
    end
    return -1
end
tspan = (0., .1)
n = parse(Int,ARGS[1])

lr= parse(Float64,ARGS[3])
AD = Optimization.AutoForwardDiff()
if ARGS[4] in ["Reverse","Rev","R","rev","r","reverse"]
    AD = Optimization.AutoReverseDiff()
end
integrator = eval(Meta.parse(ARGS[5]))
println(AD)
println(integrator)

println("$(n)mer")
nmer = get_fc_rn(n)


params = fill(10.0,n-1)
monomer_conc = fill(100.0,n)
u0 = get_species_conc(monomer_conc,nmer)

println("Iteration  t95")

for iters in 1:1000
    new_params = optim(nmer,tspan,params,monomer_conc,lr,iters,AD,integrator)
    ode = ODEProblem(nmer, u0, (.00000001,10000), new_params)
    sol = solve(ode,integrator)
    

    println(iters, " ",t95(sol,new_params))
    flush(stdout)
end
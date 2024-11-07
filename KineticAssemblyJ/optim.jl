using Catalyst
using DifferentialEquations
using Optimization
using OptimizationOptimisers
using ForwardDiff
using Zygote
#using Flux
using SciMLSensitivity
using OptimizationOptimJL

#Catalyst,DifferentialEquations,Optimization,OptimizationOptimisers,ForwardDiff,Zygote,Flux,SciMLSensitivity,OptimizationOptimJL
#for now just going to plug in max yield as lowest count monomer
function optim(rn,tspan,p_init,u0) begin
    prob = ODEProblem(rn, u0, tspan, p_init)
    p_init_values = [p[2] for p in p_init]
    max_yield = 100. # can change this later - this should determine max possible yield from the concentrations of all species -> can probably look a python code to see how this is calculated
    lr = .01 # probably make this an input parameter
    #get total yield from end of simulation for loss
    function loss(p::AbstractVector{T}, _) where T
        #Get params from rn struc -> I think p is giving me updated params automatically
        #need to get updated params for each iteration
        #David calculated koff from kon, but change C0, which shouldnt be done
        
        #Remaking ode with updated parameters
        newprob = remake(prob; p=p)

        #need to save at least at beginning and final time point to get final yield
        sol = Array(solve(newprob, Rodas5P(); saveat=tspan, maxiters=1e7))
        #what does sol give me -> want to get final yields from it
        #to get final_yield [ABC]final/max_possible[ABC]
        final = sol.u[end][end]

        yield = final/max_yield
        forward_rates = [rx.rate for rx in reactions(rn)]

        #this penalty should work for homo and hetero rates
        penalty = sum(relu.((10*lr)-forward_rates)) + sum(relu.(forward_rates-10)) 


        loss = yield + penalty

        return loss,sol

    end
    
    optf = OptimizationFunction(loss, Optimization.AutoForwardDiff())
    optprob = OptimizationProblem(optf, p_init_values)
    

    sol = solve(optprob, OptimizationOptimisers.Adam(lr); maxiters=1000)

    return sol.u
end
end
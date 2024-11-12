using Catalyst
using DifferentialEquations
using Optimization
using OptimizationOptimisers
using ForwardDiff
using Zygote
using Flux
using SciMLSensitivity
using OptimizationOptimJL

# TODO: Try using lower lr and using like 10-100 iters
# TODO: Try using forwarddiff and reversediff
# TODO: Try calculating backward rates from forward rates - can just input only forward params and calculate backwards params each iteration

#Catalyst,DifferentialEquations,Optimization,OptimizationOptimisers,ForwardDiff,Zygote,Flux,SciMLSensitivity,OptimizationOptimJL
#for now just going to plug in max yield as lowest count monomer
function optim(rn,tspan,p_init,u0) begin
    init_rates = get_rates(p_init)
    print(init_rates)
    prob = ODEProblem(rn, u0, tspan, init_rates)
    # super hard coded the max yield - just gonna assume only monomers at start and can use minimum conc of monomers
    max_yield = 100. # can change this later - this should determine max possible yield from the concentrations of all species -> can probably look a python code to see how this is calculated
    lr = .01 # probably make this an input parameter

    #get total yield from end of simulation for loss
    function loss(p::AbstractVector{T}, _) where T
        #Get params from rn struc -> I think p is giving me updated params automatically
        #need to get updated params for each iteration
        #David calculated koff from kon, but change C0, which shouldnt be done
        
        #so p is a vector of duals and so need to access value of each dual to get current parameters
        #rates = get_rates([i.value for i in p],n)
        rates = get_rates(p)

        #Remaking ode with updated parameters
        #  see how rates being float vs dual affects results - tried and doesnt make a difference
        newprob = remake(prob; p=rates)

        #need to save at least at beginning and final time point to get final yield
        sol = solve(newprob, Rodas5P(); saveat=tspan, abstol=1e-10, reltol=1e-8, maxiters=1e7)
        #sol = Array(solve(newprob, TRBDF2(); saveat=tspan, maxiters=1e7))
        #what does sol give me -> want to get final yields from it
        #to get final_yield [ABC]final/max_possible[ABC]
        final = sol.u[end][end]
        #Final value is very low for first iteration for some reason???
        #println("Final")
        #println(final)
        yield = final#/max_yield
        
        # println("Max yield")
        # println(max_yield)
        # println("Yield")
        # println(yield)
        #println("Iter: ",cnt)
        
        #this penalty should work for homo and hetero rates
        #if rate gets too high then have to take into account diffusion - makes problem PDE not ODE
        #higher on rate means low energetic barrier so problem becomes more diffusion based
        penalty = sum(relu.((10*lr).-p)) + sum(relu.(p.-10)) 

        #I am not sure if here what sign the loss should have, since using an optfxn, but David seemed to just pass normal MSE
        # println("Penalty")
        # println(penalty)
        loss = -yield + penalty
        # println("Loss")
        # println(loss)
        return loss # I think issue is type I am returning here

    end
    
    
    optf = OptimizationFunction(loss, Optimization.AutoForwardDiff())
    optprob = OptimizationProblem(optf, p_init)
    
    # TODO: Use various built in solvers 
    # actually I dont think this does integration, just Optimization
    # This does not work if I give it an integration algorithm - so not taking derivative
    # Can try inputting a integration method to see if it uses it
    sol = solve(optprob, OptimizationOptimisers.Adam(lr); maxiters=5000)

    return get_rates(sol.u)
end
end

#For fully connected m should be amount of monomers in each species
#I think C0 should be a different value maybe
function get_rates(forward_rates::Vector{T}; delta_G_kb_T::Float64=-20., C0::Float64=1e6) where T
    rates = Vector{T}()
    for (m,k_on) in enumerate(forward_rates)
        push!(rates, k_on)
        push!(rates, k_on * C0 * exp(m * delta_G_kb_T))
    end
    
    return rates
end
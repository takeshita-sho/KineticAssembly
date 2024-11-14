using Catalyst
using DifferentialEquations
using Optimization
using OptimizationOptimisers
using ForwardDiff
using Zygote
using Flux
using SciMLSensitivity
using OptimizationOptimJL
using Combinatorics: combinations

# TODO: Try using lower lr
# TODO: Try using forwarddiff and reversediff

function optim(rn,tspan,p_init,monomer_conc,lr,iters) begin
    k_symbols = unique!(reactionrates(rn))
    
    init_rates = get_rates(p_init,k_symbols)
    u0 = get_species_conc(monomer_conc)
    
    println(init_rates)
    #println("Loss     Yield")

    prob = ODEProblem(rn, u0, tspan, init_rates)

    #get total yield from end of simulation for loss
    function loss(p::AbstractVector{T}, _) where T
        #Get params from rn struc -> p gives updated params

        #David calculated koff from kon, but changes C0, which shouldnt be done
        
        #p is a vector of duals
        
        rates = get_rates(p,k_symbols)
        
        #Remaking ode with updated parameters
        newprob = remake(prob; p=rates)

        #need to save at least at beginning and final time point to get final yield
        sol = solve(newprob, Rodas5P(); saveat=[tspan[2]], abstol=1e-10, reltol=1e-8, maxiters=1e7)
        #sol = Array(solve(newprob, TRBDF2(); saveat=tspan, maxiters=1e7))
        
        #to get final_yield [ABC]final/max_possible[ABC]
        yield = sol.u[end][end]
        
        
        #this penalty should work for homo and hetero rates
        #if rate gets too high then have to take into account diffusion - makes problem PDE not ODE
        #higher on rate means low energetic barrier so problem becomes more diffusion based
        penalty = sum(relu.((10*lr).-p)) + sum(relu.(p.-10)) 

        #David seemed to just pass normal MSE
        
        loss = -yield + penalty
        
        #println(loss.value, " " ,yield.value)
        return loss

    end
    
    #optf = OptimizationFunction(loss, Optimization.AutoZygote())
    optf = OptimizationFunction(loss, Optimization.AutoForwardDiff())
    optprob = OptimizationProblem(optf, p_init)
    
    # TODO: Use various built in optimizers 
    # This does not do integration, just Optimization
    # This does not work if I give it an integration algorithm - so not taking derivative
    
    #sol = solve(optprob, Optimization.LBFGS(); maxiters=10000)
    sol = solve(optprob, OptimizationOptimisers.Adam(lr); maxiters=iters)

    return get_rates(sol.u, k_symbols)
end
end


#For fully connected m should be amount of monomers in each species
#I think C0 should be a different value maybe
function get_rates(forward_rates::Vector{T},k_symbols; delta_G_kb_T::Float64=-20., C0::Float64=1e6) where T
    rates = Vector{T}()
    for (m,k_on) in enumerate(forward_rates)
        #push!(rates, k_on)
        #push!(rates, k_on * C0 * exp(m * delta_G_kb_T))

        # This is more inefficient, but inplace modifications are not allowed with AutoZygote
        rates = [rates; [k_on, k_on * C0 * exp(m * delta_G_kb_T)]]
    end

    return dict(zip(rates,k_symbols))
end

function get_species_conc(monomers)
    n = length(monomers)
    species = Vector{Symbol}()
    
    # Generate monomer species names
    for i in 1:n
        push!(species, Symbol("X$(i)"))
    end
    
    # Generate all possible combinations for higher-order species
    for size in 2:n
        for combo in combinations(1:n, size)
            name = join(["X$i" for i in combo])
            push!(species, Symbol(name))
        end
    end
    
    # Initialize dictionary with zeros
    conc_dict = Dict{Symbol,Float64}()
    for s in species
        conc_dict[s] = 0.0
    end
    
    # Map provided monomer concentrations
    for i in 1:n
        conc_dict[Symbol("X$i")] = monomers[i]
    end
    
    return conc_dict
end
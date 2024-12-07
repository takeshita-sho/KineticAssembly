using Catalyst
using OrdinaryDiffEq, DiffEqCallbacks
using OptimizationOptimisers
using NNlib: relu
using ReverseDiff
include("./ReactionNetwork.jl")

"""
Optimization function that optimizes the rate constants of the reaction network
given initial forward rates and monomer concentrations.
"""
function optim(rn::ReactionSystem,tspan::Tuple{Float64,Float64},p_init::Vector{Float64},
    monomer_conc::Vector{Float64},lr::Float64,iters::Int,AD::AbstractADType,integrator::Any;verbose::Bool=false)
    pmap = Catalyst.paramsmap(rn)
    k_symbols = collect(keys(pmap))
    sort!(k_symbols, by = k -> pmap[k])
    
    init_rates = get_rates(p_init,k_symbols)
    u0 = get_species_conc(monomer_conc,rn)
    if verbose
        println(init_rates)
    end
    prob = ODEProblem(rn, u0, tspan, init_rates; jac = true)
    
    optf = OptimizationFunction(loss, AD)
    optprob = OptimizationProblem(optf, p_init, (k_symbols, prob, tspan, integrator, lr))
    sol = solve(optprob, OptimizationOptimisers.Adam(lr); maxiters=iters)
   
    return get_rates(sol.u, k_symbols)
end

"""
Calculates the loss by calculating the back rates from the forward rates to integrate
the ODE to calculate yield, then applies a penalty to the loss based on the forward rates.
"""
function loss(p::AbstractVector{T}, tuple::Any) where T
    k_symbols, prob, tspan, integrator, lr = tuple
    #David calculated koff from kon, but changes C0, which shouldnt be done
    
    #Remaking ode with updated parameters
    rates = get_rates(p,k_symbols)
    
    
    newprob = remake(prob;p=rates)#deepcopy(rates))

    
    sol = solve(newprob, integrator; saveat=[tspan[2]], abstol=1e-10, reltol=1e-8, maxiters=1e7)#, dtmin=1e-12, force_dtmin=true)
    
    yield = sol.u[end][end]
    
    
    #this penalty should work for homo and hetero rates
    #if rate gets too high then have to take into account diffusion - makes problem PDE not ODE
    #higher on rate means low energetic barrier so problem becomes more diffusion based
    penalty = sum(relu.((10*lr).-p)) + sum(relu.(p.-10)) 

    #David seemed to just pass normal MSE
    
    loss = -yield + penalty
    
    #println(loss.value, " " ,yield.value)
    #println(Sys.total_memory())
    
    #peak_memory_bytes = Sys.maxrss()
    #println("$(peak_memory_bytes / (1024^2)) MB")
    #println("Loss $(loss)")
    #flush(stdout)
    
    return loss

end

"""
Calculates the backwards rates based off the forward rates and then returns a dictionary
mapping the rates to the rate constant variables
"""

#For fully connected m should be amount of monomers in each species
function get_rates(forward_rates::AbstractArray{T},k_symbols::Vector{SymbolicUtils.BasicSymbolic{Real}}; delta_G_kb_T::Float64=-20., C0::Float64=1e6) where T
    rates = Dict{Num,T}(k_symbols .=> zeros(Float64, length(k_symbols)))
    #rates = Vector{T}(undef, 2 * length(forward_rates))
    for (m,k_on) in enumerate(forward_rates)
        rates[k_symbols[2m - 1]] = k_on
        rates[k_symbols[2m]] = k_on * C0 * exp(m * delta_G_kb_T)

        # This is more inefficient, but inplace modifications are not allowed with AutoZygote
        #rates = [rates; [k_on, k_on * C0 * exp(m * delta_G_kb_T)]]
    end

    return rates
end

"""
Creates a dictionary mapping the input vector concentrations to the monomers in
the reaction network and then initializes the rest of the species to 0
"""

function get_species_conc(monomers::Vector{Float64},rn::ReactionSystem)
    t = Catalyst.DEFAULT_IV
    mon_spec = [funcsym(Symbol("X",i),t) for i in 1:length(monomers)]
    
    species = Catalyst.get_species(rn)
    conc_dict = Dict{Num,Float64}()
    for (i,spec) in enumerate(mon_spec)
        conc_dict[spec] = monomers[i]
    end
    for spec in species
        if !haskey(conc_dict,spec)
            conc_dict[spec] = 0.0
        end
    end
    
    
    return conc_dict
end
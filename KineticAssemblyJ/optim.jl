using Catalyst
using OrdinaryDiffEq, DiffEqCallbacks
using Optimization
using OptimizationOptimisers
using ForwardDiff
using Zygote
using Flux
using OrderedCollections
using SciMLSensitivity
using OptimizationOptimJL
using SciMLStructures: Tunable, replace, replace!
using SymbolicIndexingInterface: parameter_values
include("./ReactionNetwork.jl")
#parameters()

#using SymbolicIndexingInterface
#pgetter = getp(odeprob) #have to insert symbols for rates
#psetter = setp(odeprob, [:a, :b, :c])

# TODO: Try using lower lr
# TODO: Try using forwarddiff and reversediff

function optim(rn,tspan,p_init,monomer_conc,lr,iters,AD,integrator;verbose=false) begin
    pmap = Catalyst.paramsmap(rn)
    pkeys = collect(keys(pmap))
    k_symbols = sort(pkeys, by = k -> pmap[k])
    
    #I think I need to feed in the pmap along with the unorganized key array so
    # I can return rates in the same order as the keys - there will be duplicate rates
    #println("getting rates")
    #flush(stdout)
    init_rates = get_rates(p_init,k_symbols)
    #println("rates gotten")
    #println("getting species conc")
    #flush(stdout)
    u0 = get_species_conc(monomer_conc,rn)
    #println("species conc gotten")
    #flush(stdout)
    if verbose
        println(init_rates)
    end
    #println("Constructing ODE with jacobian")
    #flush(stdout)
    #println("u0: ",u0)
    #println("Equations: ",ModelingToolkit.get_eqs(rn))
    #println("Unknown: ",ModelingToolkit.get_unknowns(rn))
    prob = ODEProblem(rn, u0, tspan, init_rates; jac = true)
    #println("prob.u0", prob.u0)
    #println("prob.p", prob.p)
    #println("ODE constructed")
    #flush(stdout)
    #get total yield from end of simulation for loss
    function loss(p::AbstractVector{T}, _) where T
        #Get params from rn struc -> p gives updated params

        #David calculated koff from kon, but changes C0, which shouldnt be done
        #p is a vector of duals
        #println("Loss being calculated")
        rates = get_rates(p,k_symbols)
        #println(rates)
        #Remaking ode with updated parameters
        # newprob = remake(prob; p=rates)
        #no paramsmap, params, p_init,parameters
        #println(prob.p)
        #println(parameters(rn).p)
        #newprob = remake(prob;p=rates,u0=u0)#replace(Tunable(), prob.p, rates))#p=copy(rates),u0=copy(prob.u0))#
        newprob = remake(prob;p=rates)#replace(Tunable(), prob.p, values(rates)))
        #println("p values: ",newprob.p)
        #println("rates: ",rates)
        #replace!(Tunable(), newprob.p, rates)#newprob.ps, rates)#does this need to be a vector?

        # Set new parameters
        #psetter(newprob, rates)


        #need to save at least at beginning and final time point to get final yield
        sol = solve(newprob, integrator; saveat=[tspan[2]], abstol=1e-10, reltol=1e-8, maxiters=1e7)#, dtmin=1e-12, force_dtmin=true)
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
        #println(Sys.total_memory())
        #peak_memory_bytes = Sys.maxrss()
        #println("Peak memory usage so far: $(peak_memory_bytes / (1024^2)) MB")
        #println("Loss $(loss)")
        #flush(stdout)
        return loss

    end
    
    #println("constructing optimization function")
    #flush(stdout)
    #optf = OptimizationFunction(loss, Optimization.AutoZygote())
    optf = OptimizationFunction(loss, AD)
    #println("optf constructed")
    #println("Constructing optproblem")
    #flush(stdout)
    optprob = OptimizationProblem(optf, p_init)
    #println("optp constructed")
    #flush(stdout)
    # TODO: Use various built in optimizers 
    # This does not do integration, just Optimization
    # This does not work if I give it an integration algorithm - so not taking derivative
    
    #sol = solve(optprob, Optimization.LBFGS(); maxiters=iters)
    #println("Solving")
    #flush(stdout)
    sol = solve(optprob, OptimizationOptimisers.Adam(lr); maxiters=iters)
    #println("Solved")
    #flush(stdout)
    return get_rates(sol.u, k_symbols)
end
end


#For fully connected m should be amount of monomers in each species
#I think C0 should be a different value maybe
function get_rates(forward_rates::Vector{T},k_symbols; delta_G_kb_T::Float64=-20., C0::Float64=1e6) where T
    rates = Dict{Num,T}(k_symbols .=> zeros(Float64, length(k_symbols)))
    #rates = Vector{T}(undef, 2 * length(forward_rates))
    for (m,k_on) in enumerate(forward_rates)
        rates[k_symbols[2m - 1]] = k_on
        rates[k_symbols[2m]] = k_on * C0 * exp(m * delta_G_kb_T)

        # This is more inefficient, but inplace modifications are not allowed with AutoZygote
        #rates = [rates; [k_on, k_on * C0 * exp(m * delta_G_kb_T)]]
    end

    return rates#Dict(zip(k_symbols,rates))#rates
end

#Fix this - probably indexing needs to use funsym(symbol)
function get_species_conc(monomers,rn)
    t = Catalyst.DEFAULT_IV
    mon_spec = [funcsym(Symbol("X",i),t) for i in 1:length(monomers)]
    #n = length(monomers)
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
    
    # Map provided monomer concentrations
    
    
    return conc_dict
end
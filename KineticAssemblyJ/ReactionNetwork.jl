using Catalyst
using OrderedCollections
using Catalyst.Combinatorics

"""
creates a ModelingToolkit function-like Symbol
Credit to @isaacsas for this function
"""
function funcsym(S::Symbol, t, args...)
    S = Symbol(S,args...)
    (@species $(S)(t))[1]
end

"""
Generates rate constants for a reaction network of size n
"""

function gen_rates(n)
    ptoids = OrderedDict{Symbol,Int}()
    for i in 1:(n-1)*2
        psym = Symbol("k",i)
        ptoids[psym] = i
    end
    return [(@parameters $psym)[1] for psym in keys(ptoids)]
end

"""
Generates a reaction network for a fully connected nmer
"""
function get_fc_rn(n;t=Catalyst.DEFAULT_IV)
    rxs=[]
    cnt = 1
    rates = gen_rates(n)
    for size in 2:n
        for combo in combinations(1:n, size)
            name = join(["X$i" for i in combo])
            name = funcsym(Symbol(name),t)
            for s in combo
                s1 = funcsym(Symbol("X",s),t)
                s2 = join(["X$i" for i in setdiff(combo,[s])])
                s2 = funcsym(Symbol(s2),t)
                @parameters t
                push!(rxs,Reaction(rates[cnt],[s1,s2], [name]))
                push!(rxs,Reaction(rates[cnt+1],[name], [s1,s2]))
                if length(combo) == 2
                    break
                end
            end
            
        end
        cnt+=2
    end
    @named rn = ReactionSystem(rxs,t)
    return complete(rn)
end







#=using Catalyst
using DifferentialEquations

"""
This file contains functions needed for the construction of fully connected hetero-n-mer
reaction networks as objects of the ReactionSystem type provided by Catalyst.jl.
It makes heavy use of meta-programming.
"""


"""
    get_fc_species_list(n)

Given integer `n`, return a length-sorted Vector{Vector{String}} object in which
Vector{String} objects are all ways of listing a non-empty subset of the `n` 
monomers in a fully connected reaction network topology, up to permutation.
Each monomer is represented by a string of the form `Xi`.
"""
function get_fc_species_list(n::Int64)::Vector{Vector{String}}
    monomers = ["X$i" for i in 1:n]
    species_list = []

    for i in 1:(2 ^ n - 1)
        # `inclusion` uses binary encodings to iterate through non-empty elements
        # of the power set of the set of monomers.
        inclusion = digits(i, base=2, pad=n)
        addend_species = [monomers[i] for i in 1:n if inclusion[i] == 1]
        sort!(addend_species, by=x->parse(Int, x[2:end]))
        species_list = [species_list; [addend_species]]
    end

    sort!(species_list, by=x->length(x))
    return species_list
end

"""
    get_species_string(species_list)

Given Vector{Vector{String}} `species_list`, return a String to be evaluated
as code that lists all species as functions of time for Catalyst.jl's `@species`
macro. 
"""
function get_species_string(species_list::Vector{Vector{String}})::String
    # Convert Vector{String} objects that represent species into Strings
    species_string_list = [join(species) for species in species_list]
    # Prepend macros
    eval_string_line1 = "@variables t; "
    eval_string_line2 = "@species " * join(species_string_list, "(t) ") * "(t)"

    return eval_string_line1 * eval_string_line2
end

"""
    convert_vec_to_string(vec)

Given a Vector object `vec`, return a String representation of it. E.g., given
a vector containing "X1" and "X2", return String "[X1, X2]".
"""
function convert_vec_to_string(vec::Vector{String})::String
    return "[" * join(vec, ", ") * "]"
end

"""
    get_fc_rxs_eval_string(species_list)

Given a Vector{Vector{String}} `species_list`, return a String to be evaluated
as code that lists all reactions in a fully connected reaction network topology.
"""
function get_fc_rxs_eval_string(species_list::Vector{Vector{String}})::String
    # In the fully connected topology, the number of species is 2 ^ n - 1.
    num_monomers = Int(log2(length(species_list) + 1))
    # Remove the end product from species_list
    sort!(species_list, by=x->length(x))
    pop!(species_list)

    monomers = [monomer[1] for monomer in species_list[1:num_monomers]]

    rxs_eval_string = "rxs = ["
    # reactant is a vector of bound monomers
    for reactant in species_list
        k_index = length(reactant) * 2 - 1
        # monomer is a singleton vector containing a monomer
        for monomer in monomers
            # Don't record reactions in which either (1) a complex containing itself
            # or (2) the reaction between monomers has already been recorded, just
            # in the other order of reactants.
            if monomer in reactant ||
                (length(reactant) == 1 && Int(monomer[2]) < Int(reactant[1][2]))
                continue
            end

            reactants = [join(reactant), monomer]
            product_string = join(sort([reactant; monomer], by=x->parse(Int, x[2:end])))
            product_string = "[" * product_string * "]"
            reactants_string = convert_vec_to_string(reactants)
            # Add forward direction of the reaction
            rxs_eval_string *= "Reaction(k[$(k_index)], $(reactants_string), $(product_string)), "
            # Add reverse direction of the reaction
            rxs_eval_string *= "Reaction(k[$(k_index+1)], $(product_string), $(reactants_string)), "
        end
    end

    # Prepend macro that specifies the parameters, which we store in the Vector{Float64}
    # `k`
    num_params = (num_monomers - 1) * 2
    eval_string_line1 = "@parameters k[1:$(num_params)]; "
    # Append closing bracket to Vector{Reaction} object
    eval_string_line2 = rxs_eval_string * "]"

    return eval_string_line1 * eval_string_line2
end


"""
    get_fc_rn(n)

Given integer `n`, return a ReactionSystem object representing a fully connected
rate growth chemical reaction network. This function is the master organizer of
this file, calling every function except for `get_rate_constants_from_k_ons()`.
"""
function get_fc_rn(n::Int64)::ReactionSystem
    
    species_list = get_fc_species_list(n)

    # Evaluate meta-programmed strings describing the species and reactions involved
    # in the network.
    #println(species_list)
    species_eval_string = get_species_string(species_list)
    eval(Meta.parse(species_eval_string))
    rxs_eval_string = get_fc_rxs_eval_string(species_list)
    eval(Meta.parse(rxs_eval_string))

    @named rn = ReactionSystem(rxs, t)
    return complete(rn)
end
=#
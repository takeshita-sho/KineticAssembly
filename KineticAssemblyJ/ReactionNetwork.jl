using Catalyst
using Catalyst.Combinatorics: combinations

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

function gen_rates(n::Int)
    k_symbols = Vector{Num}(undef, 2 * (n-1))
    for i in 1:(n-1)*2
        psym = Symbol("k",i)
        #ptoids[psym] = i
        k_symbols[i] = (@parameters $psym)[1]
    end
    return k_symbols
end

"""
Generates a reaction network for a fully connected nmer
"""
function get_fc_rn(n::Int;t=Catalyst.DEFAULT_IV)
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
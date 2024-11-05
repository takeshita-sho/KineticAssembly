using KineticAssemblyJ
using KineticAssemblyJ.reaction_network: gtostr
using Random
using LinearAlgebra
using SymPy

function find_eq_eqns(rn)
    # ensure same elements reference the same symbols.
    sym_buf = Dict{Int, Symbol}()
    constraints = Dict{Set{Symbol}, Int}()
    for i in keys(rn.network.nodes)
        name = Set(rn.network.nodes[i]["struct"].nodes())
        n = Symbol(gtostr(rn.network.nodes[i]["struct"]))
        sym_buf[i] = n
        copies = rn.network.nodes[i]["copies"]
        if copies > 0
            constraints[name] = -1 * copies
        end
        for key in keys(constraints)
            # if key is contained in this node
            if length(intersect(key, name)) == length(key)
                constraints[key] += n
            end
        end
    end
    eqn_list = collect(values(constraints))
    for n in keys(rn.network.nodes)
        c = sym_buf[n]
        for r_set in rn.get_reactant_sets(n)
            r_tup = Tuple(r_set)
            data = rn.network.get_edge_data(r_tup[1], n)
            a = sym_buf[r_tup[1]]
            b = sym_buf[r_tup[2]]

            kon = data["k_on"]
            koff = data["k_off"]
            println("Off rates: ", koff)
            eqn = -a * b * kon + c * koff
            push!(eqn_list, eqn)
        end
    end
    return eqn_list, sym_buf
end

struct EquilibriumSolver
    rn::ReactionNetwork
    poly_system::Vector{Expr}
    symbols::Dict{Int, Symbol}

    function EquilibriumSolver(net::ReactionNetwork)
        poly_system, symbols = find_eq_eqns(net)
        new(net, poly_system, symbols)
    end

    function solve(self::EquilibriumSolver, depth::Int=0, init_val::Vector{Float64}=[], verifyBool::Bool=true)
        if depth > 100
            println("No acceptable solution found")
            return nothing
        end
        copies = collect(values(self.rn._initial_copies))
        if isempty(init_val)
            init_val = rand(length(copies)) .* maximum(copies)
        end
        solution = nothing
        try
            solution = nsolve(self.poly_system, collect(values(self.symbols)), init_val, prec=7, max_steps=1000000000, verify=verifyBool)
        catch e
            solution = solve(self, depth + 1)
        end

        return solution
    end
end


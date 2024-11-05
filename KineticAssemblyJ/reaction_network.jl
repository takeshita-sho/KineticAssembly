#module reaction_network
using Graphs#maybe use Graphs instead
using Random
using Flux
using MetaGraphs
using Random


const LOOP_COOP_DEFAULT = 1

function _equal(n1, n2)
    nm = isomorphism(categorical_node_match("label", nothing))
    int_n1 = convert_node_labels_to_integers(n1, label_attribute="label")
    int_n2 = convert_node_labels_to_integers(n2, label_attribute="label")
    return is_isomorphic(int_n1, int_n2, node_match=nm)
end

function gtostr(g::MetaGraph)
    stout = ""
    for n in vertices(g)
        stout *= string(n)
    end
    stout = sort(collect(stout)) |> join
    return stout
end

mutable struct ReactionNetwork
    network::MetaGraph
    allowed_edges::Dict
    _node_count::Int
    _rxn_count::Int
    num_monomers::Int
    is_one_step::Bool
    rxn_coupling::Bool
    uid_map::Dict
    boolCreation_rxn::Bool
    creation_species::Vector
    creation_nodes::Vector
    creation_rxn_data::Dict
    titration_end_conc::Float64
    default_k_creation::Float64
    boolDestruction_rxn::Bool
    destruction_species::Vector
    destruction_nodes::Vector
    destruction_rxn_data::Dict
    default_k_destruction::Float64
    max_subunits::Int
    max_interactions::Int
    monomer_add_only::Bool
    chaperone::Bool
    homo_rates::Bool
    observables::Dict
    flux_vs_time::Dict
    seed::Union{Nothing, Int}
    _initial_copies::Dict
    parameters::Dict
    is_energy_set::Bool
    mon_rxns::Dict
    rxn_cid::Dict
    rxn_class::Dict
    mon_rxn_map::Dict
    dG_map::Dict
    default_k_on::Float64
    k_on::Int
    function ReactionNetwork(bngl_path::String, one_step::Bool, seed::Union{Nothing, Int}=nothing)
        network = MetaGraph(DiGraph())
        allowed_edges = Dict()
        _node_count = 0
        _rxn_count = 0
        num_monomers = 0
        rxn_coupling = false
        uid_map = Dict()
        boolCreation_rxn = false
        creation_species = Vector{Any}()
        creation_nodes = Vector{Any}()
        creation_rxn_data = Dict()
        titration_end_conc = -1.0
        default_k_creation = 1e-1
        boolDestruction_rxn = false
        destruction_species = Vector{Any}()
        destruction_nodes = Vector{Any}()
        destruction_rxn_data = Dict()
        default_k_destruction = 1e-1
        max_subunits = -1
        max_interactions = 2
        monomer_add_only = false
        chaperone = false
        homo_rates = false
        observables = Dict()
        flux_vs_time = Dict()
        _initial_copies = Dict()
        parameters = Dict()
        is_energy_set = true
        mon_rxns = Dict()
        rxn_cid = Dict()
        rxn_class = Dict()
        mon_rxn_map = Dict()
        dG_map = Dict()
    
        # Assuming parse_bngl is defined elsewhere
        #parse_bngl(open(bngl_path, "r"), seed=seed)
    
        return new(network, allowed_edges, _node_count, _rxn_count, num_monomers, 
        one_step, rxn_coupling, uid_map, boolCreation_rxn, creation_species, creation_nodes, 
        creation_rxn_data, titration_end_conc, default_k_creation, boolDestruction_rxn, 
        destruction_species, destruction_nodes, destruction_rxn_data, default_k_destruction, 
        max_subunits, max_interactions, monomer_add_only, chaperone, homo_rates, observables, 
        flux_vs_time, seed, _initial_copies, parameters, is_energy_set, mon_rxns, 
        rxn_cid, rxn_class, mon_rxn_map, dG_map)
    end
end


#=
function get_params(self::ReactionNetwork)
    keys = Vector{Any}()
    for key in self.parameters
        push!(keys,self.parameters[key])
    end
    return keys
end

function get_reactant_sets(self::ReactionNetwork, node_id::Int)
        all_predecessors = Set(self.network.in_edges(node_id))
        predecessors_vec = Vector{Any}()
        while length(all_predecessors) > 0
            found = false
            reactant = all_predecessors.pop()
            predecessors = Dict(reactant[1])
            reactant_data = self.network[reactant[1]][reactant[2]]
            poss_coreactant = nothing
            for poss_coreactant in all_predecessors
                poss_coreactant_data = self.network[poss_coreactant[1]][poss_coreactant[2]]
                if reactant_data["uid"] == poss_coreactant_data["uid"]
                    found = true
                    break
                end
            end
            if found
                delete(all_predecessors,poss_coreactant)
                push!(predecessors,poss_coreactant[1])
            end
            push!(predecessors_vec,predecessors) 
        end
        return predecessors_vec
end

function parse_param(self::ReactionNetwork, line::String)
    items = Vector{Any}(split(line,r"\s+", limit=2))
    items[2] = Meta.parse(items[2])
    print(items)
    if items[1] == "default_assoc"
        self.default_k_on = items[2]
    elseif items[1] == "rxn_coupling"
        self.rxn_coupling = items[2]
        print(self.rxn_coupling)
    elseif items[1] =="creation_rate"
        self.default_k_creation = items[2]
    elseif items[1] =="destruction_rate"
        self.default_k_destruction = items[2]
    elseif items[1] == "max_subunits"
        self.max_subunits = items[2]
    elseif items[1] == "max_interactionsa"
        self.max_interactions = items[2]
    elseif items[1] == "monomer_add_only"
        self.monomer_add_only=items[2]
    elseif items[1] == "chaperone"
        self.chaperone=items[2]
        self.chaperone_rxns = []
        self.chap_uid_map = Dict()
        self.chap_int_spec_map = Dict()
        self.optimize_species=Dict("substrate"=>[],"enz-subs"=>[])
    elseif items[1]== "homo_rates"
        self.homo_rates=items[2]
    elseif items[1]=="titration_time_int"
        print("Setting Titration End Point")
        self.titration_end_conc=items[2]
    end
    return items
end

function parse_species(self::ReactionNetwork, line::String, params::Dict)
    
    items = split(line)
    sp_info = split(items[1], r"[\),\()]")
    local init_pop
    try
        init_pop = parse(Int, items[2])
    catch e
        try
            init_pop = parse(Float64, items[2])
        catch e
            init_pop = parse(Int, params[items[2]])
        end
    end
    state_net = MetaGraph(Graph())
    set_prop!(state_net, 1, :label, sp_info[1])
    # if self.max_subunits > 0
    #     state_net = MetaGraph(Graph())
    # else
    #     state_net = Graph()
    # end
    #setindex!(state_net,sp_info[1])
    #add_vertex!(state_net, sp_info[1])
    #state_net[:label]=sp_info[1]
    set_prop!(self.network,self._node_count,:tuple,(struc=state_net, copies=[convert(Float64,init_pop)], subunits=1))
    #self.network[self._node_count] = (struc=state_net, copies=[convert(Float64,init_pop)], subunits=1)
    self._initial_copies[self._node_count] = [convert(Float64,init_pop)]
    self._node_count += 1
end




function parse_rule(self, line, params, seed=nothing, percent_negative=0.5, score_range=100)
    items = split(line, r" |, ")
    split_01 = split(items[1], "<->")
    
    if occursin("!", split_01[1])
        r_info = split(split_01[1], r"\+")
        react_1 = join(split(r_info[1], r"\\(.\!.\\)|\."))
        react_2 = join(split(r_info[2], r"|\\(.\\)|\+"))
    else
        if "null" in split_01
            # Do nothing
        else
            r_info = split(split_01[1], r"\(.+\)+\.|\(.+\)")
            print(r_info)
            react_1 = r_info[1]
            react_2 = r_info[2]
        end
    end

    self.k_on = get(params, :default_assoc, 1)
    k_off = nothing
    
    if occursin("G=", items[end])
        score = [parse(Float64, split(items[end], "=")[2])]
    else
        if seed !== nothing
            Random.seed!(seed)
        end
        score = (rand(Float64) - percent_negative) * score_range
    end
    
    if split_01[1] == "null"
        println("Found Creation rxn")
        species = split(split_01[2], r"\\(.\\)")[1]
        self.allowed_edges[Tuple(["null", species])] = [nothing, nothing, LOOP_COOP_DEFAULT, score]
        self.boolCreation_rxn = true
        push!(self.creation_species, species)
    elseif split_01[2] == "null"
        println("Found Destruction rxn")
        species = split(split_01[1], r"\\(.\\)")[1]
        self.allowed_edges[Tuple([species, "null"])] = [nothing, nothing, LOOP_COOP_DEFAULT, score]
        self.boolDestruction_rxn = true
        push!(self.destruction_species, species)
    else
        self.allowed_edges[Tuple(sort([react_1, react_2]))] = [nothing, nothing, LOOP_COOP_DEFAULT, score]
    end

    if haskey(params, :rxn_coupling)
        self.rxn_coupling = params[:rxn_coupling]
    end
end
=#

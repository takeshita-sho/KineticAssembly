import re
import sys
from typing import Tuple

import networkx as nx
import random

import torch
from torch import FloatTensor as Tensor
from torch import rand
from torch import nn

LOOP_COOP_DEFAULT = 1


def _equal(n1, n2) -> bool:
    nm = nx.algorithms.isomorphism.categorical_node_match("label", None)
    int_n1 = nx.convert_node_labels_to_integers(n1, label_attribute="label")
    int_n2 = nx.convert_node_labels_to_integers(n2, label_attribute="label")
    return nx.is_isomorphic(int_n1, int_n2, node_match=nm)
def gtostr(g: nx.DiGraph) -> str:
    stout = ""
    for n in g.nodes():
        stout += str(n)
    # make invariant
    stout = ''.join(sorted(stout))
    return stout
class ReactionNetwork:
    def __init__(self, bngl_path: str, one_step: bool, seed=None):
        self.network: nx.DiGraph() = nx.DiGraph()
        self.allowed_edges = {}
        self._node_count = 0
        self._rxn_count = 0
        self.num_monomers = 0
        self.is_one_step = one_step
        self.rxn_coupling = False
        self.uid_map = dict()
        self.boolCreation_rxn = False
        self.creation_species = []
        self.creation_nodes = []
        self.creation_rxn_data ={}
        self.titration_end_conc=-1
        self.default_k_creation = 1e-1
        self.boolDestruction_rxn = False
        self.destruction_species = []
        self.destruction_nodes = []
        self.destruction_rxn_data ={}
        self.default_k_destruction = 1e-1
        self.max_subunits = -1
        self.max_interactions = 2
        self.monomer_add_only = False
        self.chaperone=False
        self.homo_rates=False
        self.observables = dict()
        self.flux_vs_time = dict()
        self.seed = seed
        self._initial_copies = {}
        self.parse_bngl(open(bngl_path, 'r'), seed=self.seed)
        self.parameters = {} 
        self.is_energy_set = True

        self.mon_rxns = dict()
        self.rxn_cid = dict()
        self.rxn_class = dict() 
        self.mon_rxn_map = dict()
        self.dG_map = dict()

    def get_params(self):
        for key in self.parameters:
            yield self.parameters[key]

    def get_reactant_sets(self, node_id: int):
        all_predecessors = set(self.network.in_edges(node_id))
        while len(all_predecessors) > 0:
            found = False
            reactant = all_predecessors.pop()
            predecessors = {reactant[0]}
            reactant_data = self.network[reactant[0]][reactant[1]]
            poss_coreactant = None
            for poss_coreactant in all_predecessors:
                poss_coreactant_data = self.network[poss_coreactant[0]][poss_coreactant[1]]
                if reactant_data['uid'] == poss_coreactant_data['uid']:
                    found = True
                    break
            if found:
                all_predecessors.remove(poss_coreactant)
                predecessors.add(poss_coreactant[0])
            yield predecessors

    def parse_param(self, line):
        items = line.split(None, 1)
        items[1] = eval(items[1])
        print(items)
        if items[0] == 'default_assoc':
            self.default_k_on = items[1]
        elif items[0] == 'rxn_coupling':
            self.rxn_coupling = items[1]
            print(self.rxn_coupling)
        elif items[0] =='creation_rate':
            self.default_k_creation = items[1]
        elif items[0] =='destruction_rate':
            self.default_k_destruction = items[1]
        elif items[0] == 'max_subunits':
            self.max_subunits = items[1]
        elif items[0] == 'max_interactions':
            self.max_interactions = items[1]
        elif items[0] == 'monomer_add_only':
            self.monomer_add_only=items[1]
        elif items[0] == 'chaperone':
            self.chaperone=items[1]
            self.chaperone_rxns = []
            self.chap_uid_map = {}
            self.chap_int_spec_map = {}
            self.optimize_species={'substrate':[],'enz-subs':[]}
        elif items[0]== 'homo_rates':
            self.homo_rates=items[1]
        elif items[0]=='titration_time_int':
            print("Setting Titration End Point")
            self.titration_end_conc=items[1]
        return items

    def parse_species(self, line, params):
        items = line.split()
        sp_info = re.split('\\)|,|\\(', items[0])
        try:
            init_pop = int(items[1])
        except ValueError:
            try:
                init_pop = float(items[1])
            except ValueError:
                init_pop = int(params[items[1]])
        if self.max_subunits>0:
            state_net = nx.MultiGraph()
        else:
            state_net = nx.Graph()
        state_net.add_node(sp_info[0])
        self.network.add_node(self._node_count, struct=state_net, copies=Tensor([float(init_pop)]),subunits=1)
        self._initial_copies[self._node_count] = Tensor([float(init_pop)])
        self._node_count += 1

    def parse_rule(self, line, params, seed=None, percent_negative=.5, score_range=100):
        items = re.split(r' |, ', line)
        split_01 = re.split('<->',items[0])
        if '!' in split_01[0]:
            r_info = re.split('\+',split_01[0])

            react_1 = "".join(re.split('\\(.\!.\\)|\.',r_info[0]))
            react_2 = "".join(re.split('|\\(.\\)|\+',r_info[1]))

        else:
            if 'null' in split_01:
                pass
            else:
                r_info = re.split('\\(.\\)+.|\\(.\\)',split_01[0])
                react_1 = r_info[0]
                react_2 = r_info[1]



        if params['default_assoc']:
            self.k_on = params['default_assoc']
        else:
            self.k_on = 1
        k_off = None
        if 'G=' in items[-1]:
            score = Tensor([float(items[-1].split('=')[1])])
        else:
            if seed:
                torch.random.manual_seed(seed)
            score = (rand(1, dtype=torch.float) - percent_negative) * score_range

        if split_01[0]=='null':
            print("Found Creation rxn")
            species = re.split('\\(.\\)',split_01[1])[0]
            self.allowed_edges[tuple(['null',species])] = [None, None, LOOP_COOP_DEFAULT, score]
            self.boolCreation_rxn=True
            self.creation_species.append(species)
        elif split_01[1]=='null':
            print("Found Destruction rxn")
            species=re.split('\\(.\\)',split_01[0])[0]
            self.allowed_edges[tuple([species,'null'])] = [None, None, LOOP_COOP_DEFAULT, score]
            self.boolDestruction_rxn=True
            self.destruction_species.append(species)
        else:
            self.allowed_edges[tuple(sorted([react_1, react_2]))] = [None, None, LOOP_COOP_DEFAULT, score]


        if params.get('rxn_coupling') is not None:
            self.rxn_coupling=params['rxn_coupling']

    def parse_bngl(self, f, seed=None):
        parameters = dict()
        cur_block = ''
        for line in f:
            line = line.strip()
            if len(line) > 0 and line[0] != '#':
                if "begin parameters" in line:
                    cur_block = 'param'
                elif "begin species" in line:
                    cur_block = 'species'
                elif "begin rules" in line:
                    cur_block = 'rules'
                elif "begin observables" in line:
                    cur_block = 'observables'
                elif "end" in line:
                    cur_block = ' '
                else:
                    if cur_block == 'param':
                        items = self.parse_param(line)
                        parameters[items[0]] = items[1]
                    elif cur_block == 'species':
                        self.parse_species(line, parameters)
                    elif cur_block == 'rules':
                        self.parse_rule(line, parameters, seed=None)


        if "loop_coop" in parameters:
            if len(parameters['loop_coop']) != len(self.allowed_edges):
                raise ValueError('num loop_coop must equal to num allowed_edges')
            keys = list(self.allowed_edges.keys())
            for i, lcf in enumerate(parameters['loop_coop']):
                if lcf > 1 or lcf < 0:
                    raise ValueError('loop cooperativity factor must be between 0 and 1')
                self.allowed_edges[keys[i]][2] = lcf
        self.num_monomers = self._node_count

    def reset(self):
        for key in self._initial_copies:
            self.network.nodes[key]['copies'] = self._initial_copies[key]
        self.observables = {}
        self.flux_vs_time = {}
        for i in range(self.num_monomers):
            self.observables[i] = (gtostr(self.network.nodes[i]['struct']), [])
            self.flux_vs_time[i] = (gtostr(self.network.nodes[i]['struct']), [])
        fin_dex = len(self.network.nodes) - 1
        self.observables[fin_dex] = (gtostr(self.network.nodes[fin_dex]['struct']), [])
        self.flux_vs_time[fin_dex] = (gtostr(self.network.nodes[fin_dex]['struct']), [])

    def intialize_activations(self, mode="middle"):
        if not self.is_energy_set:
            raise ValueError("The network free energies must be calculated for activation params to be used")
        for node in self.network.nodes:
            for reactant_set in self.get_reactant_sets(node):
                if mode == 'uniform':
                    k_on = nn.Parameter(rand(1, dtype=torch.float) * Tensor([1]), requires_grad=True)
                elif mode == 'middle':
                    k_on = nn.Parameter(Tensor([1.]), requires_grad=True)
                self.parameters[tuple(list(reactant_set) + [node])] = k_on
                for source in reactant_set:
                    self.network.edges[(source, node)]['k_on'] = k_on

    def initialize_random_pairwise_energy(self, percent_negative=.5, score_range=1000, seed=None):
        for node in self.network.nodes:
            for reactant_set in self.get_reactant_sets(node):
                if seed is not None:
                    torch.random.manual_seed(seed)
                score = (rand(1, dtype=torch.float) - percent_negative) * score_range
                for source in reactant_set:
                    if source < self.num_monomers:
                        self.network.edges[(source, node)]['rxn_score'] = score

        self.is_energy_set = True

    def _add_graph_state(self, connected_item: nx.Graph, source_1: int, source_2: int = None, template=None,subunits=1):
        if type(source_1) is not int:
            source_1 = int(source_1[0])
        if source_2 is not None and type(source_2) is not int:
            source_2 = int(source_2[0])
        node_exists = [x for x in self.network.nodes(data=True) if
                       _equal(x[1]['struct'], connected_item)]
        new_edges_added = 0
        if len(node_exists) == 0:
            print("New node added - Node index: %d ; Node label: %s " %(self._node_count,gtostr(connected_item)))
            self.network.add_node(self._node_count, struct=connected_item, copies=Tensor([0.]),subunits=subunits)
            self._initial_copies[self._node_count] = Tensor([0.])
            new_node = self._node_count
            self._node_count += 1
        elif len(node_exists) > 1:
            raise Exception("Duplicate nodes in reaction Network")
        else:
            new_node = node_exists[0][0]
        if self.network.has_edge(source_1, new_node):
            return None
        if not template:
            return None
        else:

            dg_coop = sum([self.allowed_edges[tuple(sorted(e))][3] for e in template])
            self.network.add_edge(source_1, new_node,
                                  k_on=self.default_k_on,
                                  k_off=None,
                                  lcf=1,
                                  rxn_score=dg_coop,
                                  uid=self._rxn_count)
            new_edges_added+=1
            if source_2 is not None:
                self.network.add_edge(source_2, new_node,
                                      k_on=self.default_k_on,
                                      k_off=None,
                                      lcf=1,
                                      rxn_score=dg_coop,
                                      uid=self._rxn_count)
                new_edges_added+=1
                if self.rxn_coupling and (self.network.nodes[source_1]['struct'].number_of_edges() == 0 and self.network.nodes[source_2]['struct'].number_of_edges() ==0):
                    reactants = tuple(sorted((source_1,source_2))) 
                    self.mon_rxns[reactants]=self._rxn_count
                self.uid_map[self._rxn_count] = tuple(sorted((source_1,source_2)))

            if len(template) > new_edges_added:
                print("The number of bonds formed are not compensated by the number of edges")
                print("This could be possible due to presence of a repeating subunit")
                print("SOurce1: ",source_1,source_2)

                cmn_reactant = set(template[0])
                for b in range(len(template)-1):
                    cmn_reactant = cmn_reactant.intersection(template[b+1])

                if cmn_reactant:
                    cmn_reactant = cmn_reactant.pop()
                    print("The common reactant is: ",cmn_reactant)
                    cmn_node=-1
                    for node in self.network.nodes(data=True):
                        if cmn_reactant == gtostr(node[1]['struct']) and (node[1]['struct'].number_of_edges()==0):
                            
                            cmn_node = node[0]
                    print("Edge added between: ", cmn_node,new_node)
                    self.network.add_edge(cmn_node, new_node,
                                  k_on=self.default_k_on,
                                  k_off=None,
                                  lcf=1,
                                  rxn_score=dg_coop,
                                  uid=self._rxn_count)
                    new_edges_added+=1

        self._rxn_count += 1
        if len(node_exists) == 0:
            return (new_node, self.network.nodes[new_node])
        else:
            return None

    def match_maker(self, n1, n2=None, one_step=False) -> list:
        
        nodes_added = []
        orig = n1[1]['struct']

        if n2 is not None:
            nextn = n2[1]['struct']
            item = nx.compose(orig, nextn)

        else:
            item = orig
        connected_item = item.copy()
        new_bonds = []
        add_to_graph = False
        complex_size = 0
        for poss_edge in list(self.allowed_edges.keys()):
            if False not in [item.has_node(n) for n in poss_edge] and \
                    (n2 is None or
                     (True in [orig.has_node(n) for n in poss_edge] and
                      True in [nextn.has_node(n) for n in poss_edge]))\
                    and not item.has_edge(poss_edge[0], poss_edge[1]):
                repeat_units=False

                if self.monomer_add_only==True:
                    if (orig.number_of_edges() ==0 or nextn.number_of_edges() ==0):

                        if self.chaperone and True in [item.has_node(sp) for sp in list(self.chap_int_spec_map.keys())]:
                            
                            continue
                        connected_item.add_edge(poss_edge[0], poss_edge[1])
                        new_bonds.append(poss_edge)
                        complex_size += n1[1]['subunits']
                        if n2 is not None:
                            complex_size+=n2[1]['subunits']


                        
                        for (u,v) in item.edges:
                            if u==v:
                                # print("Repeat Units")
                                repeat_units = True
                        if repeat_units:
                            connected_item.add_edge(poss_edge[1], poss_edge[0])
                            new_bonds.append(poss_edge)

                        add_to_graph=True
                        
                elif self.monomer_add_only==False:
                    if self.chaperone and True in [item.has_node(sp) for sp in list(self.chap_int_spec_map.keys())]:
                        
                        continue
                    connected_item.add_edge(poss_edge[0], poss_edge[1])
                    new_bonds.append(poss_edge)
                    complex_size += n1[1]['subunits']
                    if n2 is not None:
                        complex_size+=n2[1]['subunits']


                    
                    for (u,v) in item.edges:
                        if u==v:
                            
                            repeat_units = True
                    if repeat_units:
                        connected_item.add_edge(poss_edge[1], poss_edge[0])
                        new_bonds.append(poss_edge)

                    add_to_graph=True
                else:
                    if (orig.number_of_edges() > 0 and nextn.number_of_edges() >0):
                        if self.chaperone and True in [item.has_node(sp) for sp in list(self.chap_int_spec_map.keys())]:
                            
                            continue
                        connected_item.add_edge(poss_edge[0], poss_edge[1])
                        new_bonds.append(poss_edge)
                        complex_size += n1[1]['subunits']
                        if n2 is not None:
                            complex_size+=n2[1]['subunits']

                        for (u,v) in item.edges:
                            if u==v:
                            
                                repeat_units = True
                        if repeat_units:
                            connected_item.add_edge(poss_edge[1], poss_edge[0])
                            new_bonds.append(poss_edge)

                        add_to_graph=True


            elif True in [item.has_node(n) for n in poss_edge] and (n2 is None) and item.has_edge(poss_edge[0], poss_edge[1]):
                
                new_bonds.append(poss_edge)
                complex_size+=n1[1]['subunits']   

            elif (n2 is not None) and (True in [orig.has_node(n) for n in poss_edge] and True in [nextn.has_node(n) for n in poss_edge]) and item.has_edge(poss_edge[0], poss_edge[1]):
                
                if orig.number_of_edges() ==0 or nextn.number_of_edges() ==0:
                    

                    n_edges = orig.number_of_edges() if orig.number_of_edges() else nextn.number_of_edges()
                    print(n_edges)
                    print(connected_item.edges())
                    total_subunits = n1[1]['subunits'] + n2[1]['subunits']


                    if total_subunits <= self.max_subunits:

                        complex_size=total_subunits


                        e1 = orig.number_of_edges()
                        e2 = nextn.number_of_edges()
                        reactant_set = tuple([r1 for r1 in orig.nodes()] + [r2 for r2 in nextn.nodes()])
                        if reactant_set == poss_edge:
                            
                            if self.max_interactions ==2:
                                
                                new_bonds.append(poss_edge)
                                connected_item.add_edge(poss_edge[1], poss_edge[0])
                                if total_subunits == self.max_subunits:
                                    
                                    new_bonds.append(poss_edge)
                                    connected_item.add_edge(poss_edge[1], poss_edge[0])
                                print(connected_item.edges())
                            else:
                                print("Forming bonds to achieve max interactions from each sub-unit")
                                for i in range(n1[1]['subunits']):
                                    for j in range(n2[1]['subunits']):
                                        # print("NEW BOND ADDEDDDDDDD")
                                        new_bonds.append(poss_edge)
                                        connected_item.add_edge(poss_edge[1], poss_edge[0])

                                
                        add_to_graph=True

                else:
                    
                    if self.monomer_add_only == -1:
                        
                        for edge2 in nextn.edges():
                            connected_item.add_edge(edge2[0],edge2[1])

                        total_subunits = n1[1]['subunits'] + n2[1]['subunits']
                        if total_subunits <= self.max_subunits:
                            complex_size=total_subunits
                            reactant_set = tuple([r1 for r1 in orig.nodes()] + [r2 for r2 in nextn.nodes()])
                            if reactant_set == poss_edge:
                                
                                if self.max_interactions ==2:

                                    new_bonds.append(poss_edge)
                                    connected_item.add_edge(poss_edge[1], poss_edge[0])
                                    if total_subunits == self.max_subunits:
                                        

                                        new_bonds.append(poss_edge)
                                        connected_item.add_edge(poss_edge[1], poss_edge[0])
                                    print(connected_item.edges())
                                else:

                                    for i in range(n1[1]['subunits']):
                                        for j in range(n2[1]['subunits']):

                                            new_bonds.append(poss_edge)
                                            connected_item.add_edge(poss_edge[1], poss_edge[0])

                            add_to_graph=True

            elif (True in [item.has_node(n) for n in poss_edge]) and (True in [len(n)>1 for n in poss_edge]) and self.chaperone and n2 is not None:


                
                rxn_is_possible=False
                node_labels= (gtostr(orig),gtostr(nextn))

                if set(node_labels) == set(poss_edge):
                    print("*******Chaperone Reaction**********")
                    reactants = sorted((n1[0],n2[0]))
                    products = sorted(list(item.nodes()))
                    print(reactants,products)

                    if (reactants,products) not in self.chaperone_rxns:
                        
                        self.chaperone_rxns.append((reactants,products))

                        
                        connected_item = item.copy()
                        new_bonds.append(poss_edge)

                        sp_len = [len(e) for e in poss_edge]     
                        connected_item.add_edge(poss_edge[sp_len.index(1)][0],poss_edge[sp_len.index(1)])  

                        if poss_edge[sp_len.index(1)] not in self.chap_int_spec_map:
                            self.chap_int_spec_map[poss_edge[sp_len.index(1)]] =  [self._node_count]
                        else:
                            self.chap_int_spec_map[poss_edge[sp_len.index(1)]].append(self._node_count)

                        add_to_graph=True

                continue
        
        if one_step and add_to_graph:
            new_node = self._add_graph_state(connected_item, n1, source_2=n2, template=new_bonds,subunits=complex_size)
            if new_node is not None:
                nodes_added.append(new_node)

        return nodes_added

    def is_hindered(self, n1, n2) -> bool:
        """
        Determines if binding two species would be sterically hindered.
        :param n1: node 1 (species 1)
        :param n2: node 2 (species 2)
        :return:
        """
        node_set1 = set(n1[1]['struct'].nodes())
        node_set2 = set(n2[1]['struct'].nodes())
        
        return len(node_set1 - node_set2) < len(node_set1)

    def decompose_monomers(self,n1,monomer_set):
        if len(self.network.in_edges(n1)) == 0:
            return(True,monomer_set)
        else:
            for incoming_edge in self.network.in_edges(n1):
                flag,monomer_set = self.decompose_monomers(incoming_edge[0],monomer_set)
                if flag:
                    monomer_set.append(incoming_edge[0])
            return(False,monomer_set)

    def map_coupled_rxns(self):
        cid={}
        for uid,reactants in self.uid_map.items():
            if self.network.nodes[reactants[0]]['struct'].number_of_edges() ==0 and self.network.nodes[reactants[1]]['struct'].number_of_edges() ==0:
                
                continue
            elif self.network.nodes[reactants[0]]['struct'].number_of_edges() ==0 :
                #Reactant 1 is monomer. Reactant 2 is not
                #Get all nodes of all monomer species
                flag,monomer_set = self.decompose_monomers(reactants[1],[])
                monomer_set = list(set(monomer_set))
                for mon in monomer_set:
                    rxn_pair = tuple(sorted((reactants[0],mon)))
                    if self.mon_rxns.get(rxn_pair) is not None:
                        mon_rxn_id = self.mon_rxns[rxn_pair]
                        if uid in cid.keys() :
                            if mon_rxn_id not in cid[uid]:
                                cid[uid].append(mon_rxn_id)
                        else:
                            cid[uid] = [mon_rxn_id]
            elif self.network.nodes[reactants[1]]['struct'].number_of_edges() ==0 :
                #Reactant 2 is monomer. Reactant 1 is not
                #Get all nodes of all monomer species
                flag,monomer_set = self.decompose_monomers(reactants[0],[])
                monomer_set = list(set(monomer_set))
                for mon in monomer_set:
                    rxn_pair = tuple(sorted((reactants[1],mon)))
                    if self.mon_rxns.get(rxn_pair) is not None:
                        mon_rxn_id = self.mon_rxns[rxn_pair]
                        if uid in cid.keys() :
                            if mon_rxn_id not in cid[uid]:
                                cid[uid].append(mon_rxn_id)
                        else:
                            cid[uid] = [mon_rxn_id]
            else:
                #Both reactants are not monomers
                #Get nodes of all monomer species for each complex
                flag1,monomer_set1 = self.decompose_monomers(reactants[0],[])
                flag2,monomer_set2 = self.decompose_monomers(reactants[1],[])
                monomer_set1 = list(set(monomer_set1))
                monomer_set2 = list(set(monomer_set2))

                for m1 in monomer_set1:
                    for m2 in monomer_set2:
                        rxn_pair = tuple(sorted((m1,m2)))
                        if self.mon_rxns.get(rxn_pair) is not None:
                            mon_rxn_id = self.mon_rxns[rxn_pair]
                            if uid in cid.keys():
                                if mon_rxn_id not in cid[uid]:
                                    cid[uid].append(mon_rxn_id)
                            else:
                                cid[uid]=[mon_rxn_id]
        return(cid)

    def check_if_node_exists(self,species):
        for node in self.network.nodes():
            node_label = gtostr(self.network.nodes[node]['struct'])
            if node_label == species:
                return(True)
        return(False)


    def resolve_creation_rxn(self):

        print(self.creation_species)
        for n in self.network.nodes:
            node_lb = gtostr(self.network.nodes()[n]['struct'])
            if (node_lb in self.creation_species) and self.network.nodes[n]['struct'].number_of_edges() ==0:    #The second condition is just to check if its a monomer and not a homodimer
                self.creation_nodes.append(n)
                self.creation_rxn_data[n] = {'uid':self._rxn_count,'k_on':self.default_k_creation}
                self._rxn_count+=1
            if (node_lb in self.destruction_species) and self.network.nodes[n]['struct'].number_of_edges() ==0:
                self.destruction_nodes.append(n)
                self.destruction_rxn_data[n] = {'uid':self._rxn_count,'k_on':self.default_k_destruction}
                self._rxn_count+=1

    def resolve_chaperone_rxn(self):
        print("Resolving Chaperone Rxns::")
        print(self.chaperone_rxns)
        for chap in self.chaperone_rxns:
            reactant = chap[0]    #Which nodes are reacting. e.g. AB + X
            products=[]
            enz_sub_complx = "".join(chap[1])   #Name of enzymen subtrate complex = ABX
            chap_species = -1
            for n in self.network.nodes():
                sp_label = gtostr(self.network.nodes[n]['struct'])
                if ( sp_label in chap[1]):
                    products.append(n)
                if (n in reactant) and (sp_label in list(self.chap_int_spec_map.keys())):
                    chap_species = n
                    for int_species in self.chap_int_spec_map[sp_label]:
                        if gtostr(self.network.nodes[int_species]['struct']) == enz_sub_complx:
                            r=int_species

                if (n in reactant) and len(sp_label)>1:
                    self.optimize_species['substrate'].append(n)




            self.optimize_species['enz-subs'].append(r)

            for p in products:
                self.network.add_edge(r, p,
                              k_on=self.default_k_on,
                              k_off=None,
                              lcf=1,
                              rxn_score=torch.Tensor([float(-100)]),
                              uid=self._rxn_count)

            self.uid_map[self._rxn_count] = reactant
            if chap_species not in self.chap_uid_map:
                self.chap_uid_map[chap_species] = [self._rxn_count]
            else:
                self.chap_uid_map[chap_species].append(self._rxn_count)
            self._rxn_count+=1

            for edge in self.network.in_edges(r):
                data = self.network.get_edge_data(edge[0],edge[1])
                uid = data['uid']
                if uid not in self.chap_uid_map[chap_species]:
                    self.chap_uid_map[chap_species].append(uid)

    def create_rxn_class(self):
        uid_dict = {}
        uid_reactants = {}
        for n in self.network.nodes():
            #print(n)
            #print(rn.network.nodes()[n])
            for k,v in self.network[n].items():
                uid = v['uid']
                r1 = set(gtostr(self.network.nodes[n]['struct']))
                p = set(gtostr(self.network.nodes[k]['struct']))
                r2 = p-r1
                reactants = (r1,r2)
                uid_val = {'uid':uid,'reactants':reactants,'kon':v['k_on'],'score':v['rxn_score'],'koff':v['k_off']}
                uid_reactants[uid]=reactants
                if uid not in uid_dict.keys():
                    uid_dict[uid] = uid_val

        final_rxn_class = {}
        for key,rnts in sorted(uid_reactants.items()):
        #     print(key,"\t\t",rnts)

            l1 = len(rnts[0])
            l2 = len(rnts[1])


            if (l1,l2) in final_rxn_class.keys():
                final_rxn_class[(l1,l2)].append(key)
            elif (l2,l1) in final_rxn_class.keys():
                final_rxn_class[(l2,l1)].append(key)
            else:
                final_rxn_class[(l1,l2)] = [key]
        self.rxn_class = final_rxn_class


    def resolve_tree(self):
        """
        Build the full reaction network from whatever initial info was given
        :param is_one_step:
        :return:
        """
        new_nodes = list(self.network.nodes(data=True))
        while len(new_nodes) > 0:
            node = new_nodes.pop(0)
            for anode in list(self.network.nodes(data=True)):
                # print("Node-1 : ",node)
                # print("Node-2 : ",anode)
                if not self.is_hindered(node, anode):
                    new_nodes += self.match_maker(node, anode, self.is_one_step)
                else:
                    
                    if (node[1]['subunits']+anode[1]['subunits'] <= self.max_subunits) and (self.max_subunits >0):
                        print("Adding another subunit")
                        new_nodes+= self.match_maker(node,anode,self.is_one_step)
                    elif (node[1]['subunits']+anode[1]['subunits'] > self.max_subunits) and (self.max_subunits >0):
                        print("Max subunits limit reached")
                        # print(node[1]['struct'].edges())
                        # print(anode[1]['struct'].edges())


            # must also try internal bonds
            print("Trying internal bonds")
            new_nodes += self.match_maker(node,one_step=self.is_one_step)

        #Calculating dG of final complex
        #Add code here

        # add default observables
        #Add all nodes as observables
        for i in range(len(self.network.nodes)):
            self.observables[i] = (gtostr(self.network.nodes[i]['struct']), [])
            self.flux_vs_time[i] = (gtostr(self.network.nodes[i]['struct']), [])
        # fin_dex = len(self.network.nodes) - 1
        # self.observables[fin_dex] = (gtostr(self.network.nodes[fin_dex]['struct']), [])

        #Create rxn class dict; Used for parametrizing homogeneous model
        self.create_rxn_class()

        if self.rxn_coupling:
            self.rxn_cid = self.map_coupled_rxns()
            print("Coupling Reaction ID: ", self.rxn_cid)
        if self.boolCreation_rxn or self.boolDestruction_rxn:
            print("Resolving Creation and Destruction rxns")
            self.resolve_creation_rxn()
            print("Creation Reactions: ")
            print(self.creation_nodes)
            print(self.creation_rxn_data)
            print("Destructions Reactions: ")
            print(self.destruction_nodes)
            print(self.destruction_rxn_data)

        if self.chaperone:
            self.resolve_chaperone_rxn()

        print("Reaction Network Completed")

if __name__ == '__main__':
    bngls_path = sys.argv[1]  # path to bngl
    dt = float(sys.argv[2])  # time step in seconds
    iter = int(sys.argv[3])  # number of time steps to simulate
    m = ReactionNetwork(sys.argv[1],True)
    print('done')

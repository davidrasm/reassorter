##########################
##
## Modified from Espalier: A Python package for tree reconciliation and reconstructing ARGs using maximum agreement forests.
##
## Copyright 2021-2022 David A. Rasmussen (drasmus@ncsu.edu)
##
############################

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar
import pandas as pd
import ast

class SCAR(object):
    
    """
        Model class for the structured coalescent with ancestral recombination (SCAR)
        Ancestral states can be given (known) or marginalized over (unknown).
        This version relaxes rules about the ARG:
        Specifically coal nodes need not be strictly bifurcating and can have any number of children.
        But all children must have one and only one parent.
    """
    
    def __init__(self,reassortment_rate,M,Ne,genome_length,**kwargs):
        
        '''             
            Parameters: 
                reassortment_rate (float): reassortment rate per lineage per segment
                M (2D list/array): migration rates (forward-time) between subpopulations 
                Ne (1D list/array): effective population size of each subpopulation
                genome_length (int): The number of segments K, with unit segments [0,1), ..., [K-1, K).
               
            Optional keyword arguments:
                bounds (tuple): lower and upper bounds on estimated parameter given as (lower,upper)
                dt_step (str): integration time step used to compute lineage state probabilities
                known_ancestral_states (boolean): ancestral states need to be given in tree series (ts) if True
        '''
        
        # Model params
        self.reassortment_rate = reassortment_rate
        self.M = np.transpose(np.array(M)) # transpose to get reverse time matrix        
        self.Ne = np.array(Ne)
        self.genome_length = genome_length
        
        # Likelihood calculation params
        self.bounds = kwargs.get('bounds', (0, np.inf))
        self.dt_step = kwargs.get('dt_step', 0.1)
        self.known_ancestral_states = kwargs.get('known_ancestral_states', False)
        
    def compute_neg_log_like(self,reassortment_rate,table_nodes,table_edges):
        
        """
            Compute the negative log likelihood of ARG in TreeSequence under the SCAR model
            Here we compute the negatve log likelihood because the scipy opt minimizes this func
            Note: ts, M, Ne, genome_length need to be passed through optimizer as a tuple
            
            Parameters:
                reassortment_rate (float): reassortment rate per lineage per segment
                table_nodes (string): path to CSV containing ARG nodes
                table_edges (string): path to CSV containg ARG edges

                
        """
    
        # Hard-coding number of pops as 1
        pops = 1
        if pops == 0: # we never assigned populations
            pops = 1

        # Get transition rate matrix
        Q = self.M - np.diag(np.sum(self.M,axis=1)) # set diagonals to negative row sums

        # Read in tables of nodes and edges for ARG
        # The node id must equal its row position after sorting by time ascending (as is the default for tree sequence object)
        nodes_df = pd.read_csv(table_nodes)
        edges_df = pd.read_csv(table_edges)
        
        children = edges_df["child"].to_numpy()
        parents = edges_df["parent"].to_numpy()
        segments = edges_df["segments"].apply(ast.literal_eval).to_numpy()
 
        
        # Init lineage arrays
        active_lines = [] # active lines in ARG
        active_segments = [] # tracks ancestral segments that can reassort on active_lines
        line_state_probs = [] # lineage state probabilities for active lines
        log_like = 0.0 # log likelihood of full tree
        
        # Iterate through each event in ARG nodes working backwards through time
        for i in np.arange(0,len(nodes_df)):
            
            # Get time of event and time of next event
            event_time = nodes_df.time[i]
            if (i+1 < len(nodes_df)): # if not at final event
                next_time = nodes_df.time[i+1]
            else:
                next_time = nodes_df.time[i]
            t_elapsed = next_time - event_time # time elapsed between events

            # Determine event type from nodes_df.type 
            event_type = None
            if nodes_df.type[i] == 'sample':
                event_type = 'sample'
            if nodes_df.type[i] == 'coalescence':
                event_type = 'coalescent'
            if nodes_df.type[i] == 'recombination':
                event_type = 'recombination' # reassortment, kept as 'recombination' to match node type label
            #if event.flags == 262144:
            #    event_type = 'hidden_coalescent'
            #if event.flags == 524288:
            #    event_type = 'migration'
            
            # Initialize prob of observing events or no events
            event_prob = 1.0
            prob_no_coal = 1.0
            prob_no_mig = 1.0
            prob_no_reassort = 1.0
            
            # Update active lineages based on event type: coalescent/sampling/migration events
            if 'sample' == event_type:
                
                # Add sampled lineage
                active_lines.append(i)
                # Calculate how many ancestral segments the active lineage carries
                active_segments.append(self._get_line_segment_count(i,children,segments))
                state_probs = np.zeros(pops)
                #if event.population == -1: # we never assigned populations
                state_probs[0] = 1.0 # set prob to 1.0 for sampled state
                #else:
                #    state_probs[event.population] = 1.0 # set prob to 1.0 for sampled state
                line_state_probs.append(state_probs)            
            
            if 'coalescent' == event_type:
                
                # Get children of parent node at coalescent event
                coal_children = children[parents == i] # parent has id == i in parent column of edges table
                
                # Get uniique children b/c the same parent/child edge may occur more than once in the tree series if not in contiguous local trees
                coal_children = np.unique(coal_children)
                
                # Find coal_children in active_lines
                coal_children = [x for x in coal_children if x in active_lines]
                child_indexes = [active_lines.index(x) for x in coal_children]
                
                # Compute coalescent event prob for arbitrary number of children 
                coal_probs = np.ones(pops)
                for child_idx in child_indexes:
                    coal_probs *= line_state_probs[child_idx]
                coal_probs = coal_probs / self.Ne
                lambda_sum = sum(coal_probs)
                event_prob = lambda_sum
                
                # Compute new parent state probs
                #if self.known_ancestral_states:
                #    parent_probs = np.zeros(pops)
                #    parent_probs[event.population] = 1.0
                #else:
                parent_probs = coal_probs / lambda_sum # renormalize probs
                    
                # Update lineage arrays - overwriting child1 with parent
                active_lines[child_indexes[0]] = i # name of parent
                active_segments[child_indexes[0]] = self._get_line_segment_count(i,children,segments)
                line_state_probs[child_indexes[0]] = parent_probs
                child_indexes.pop(0) # remove first index given to parent
                for child_idx in sorted(child_indexes, reverse=True): # remove in reverse order so indexes don't change
                    del active_lines[child_idx]
                    del active_segments[child_idx]
                    del line_state_probs[child_idx]
            
            #if 'hidden_coalescent' == event_type:
                
                # Hidden coalescent in ARG not observed in local trees - only need to update active_lines but nothing else
                
            #    coal_children = children[parents == i]
            #    coal_children = np.unique(coal_children)
            #    child1 = coal_children[0]
            #    child2 = coal_children[1]
            #    child1_idx = active_lines.index(child1)
            #    child2_idx = active_lines.index(child2)
                
                # Compute likelihood of coalescent event
            #    p1 = line_state_probs[child1_idx]
            #    p2 = line_state_probs[child2_idx]
            #    coal_probs = (p1 * p2) / self.Ne
            #    lambda_sum = sum(coal_probs)
            #    event_prob = lambda_sum
                
                # Compute new parent state probs"
            #    if self.known_ancestral_states:
            #        parent_probs = np.zeros(pops)
            #        parent_probs[event.population] = 1.0
            #    else:
            #        parent_probs = coal_probs / lambda_sum
                
                # Update lineage arrays - overwriting child1 with parent"
            #    active_lines[child1_idx] = i # name of parent
            #    active_segments[child1_idx] = self._get_line_segment_count(i,children,segments)
            #    line_state_probs[child1_idx] = parent_probs
            #    del active_lines[child2_idx]
            #    del active_segments[child2_idx]
            #    del line_state_probs[child2_idx]
            
            if 'recombination' == event_type:
                
                # Find child of parent node 
                child = children[parents == i]
                child = np.unique(child)
                assert len(child) == 1
                child = child[0]

                # Remember that child may have already been removed from active_lines
                if child in active_lines:
                    
                    # Node i's two parents are what this event splits into
                    recomb_parents = parents[children == i]
                    recomb_parents = np.unique(recomb_parents)
                    
                    child_idx = active_lines.index(child)
                    
                    if len(recomb_parents) == 1:
                        # All segments went to the same parent by chance, so reassortment is unobservable
                        # This does not contribute to the likelihood
                        parent = recomb_parents[0]
                        active_lines[child_idx] = parent
                        active_segments[child_idx] = self._get_line_segment_count(parent, children, segments)
                        event_prob = 1.0
                    
                    else:
                        assert len(recomb_parents) == 2, \
                            f"Node {i}: expected 1 or 2 parents, got {len(recomb_parents)}"
                        # If there are two parents, this is a detectable reassortment event
                        # This does contribute to the likelihood 

                        left_parent, right_parent = recomb_parents[0], recomb_parents[1]
        
                        links = active_segments[child_idx]  # segment count of the pre-split lineage
                        event_prob = reassortment_rate * (1 - (0.5)**(links - 1))
                        
                        parent_probs = line_state_probs[child_idx]
                        
                        # Isolate each branch's own edge so you can track which segments correspond to a lineage
                        left_mask = (parents == left_parent)
                        right_mask = (parents == right_parent)
                        
                        # Relabel both branches with i (this node's own id) - matches the
                        # convention used by the coalescent block
                        active_lines[child_idx] = i
                        active_segments[child_idx] = self._get_line_segment_count(
                            i, children[left_mask], segments[left_mask]
                        )
                        line_state_probs[child_idx] = parent_probs
                        
                        active_lines.append(i)
                        active_segments.append(self._get_line_segment_count(
                            i, children[right_mask], segments[right_mask]
                        ))
                        line_state_probs.append(parent_probs)
            
            #if 'migration' == event_type:
                
                # Find migrating (child) lineage
            #    mig_child = children[parents == i] # parent has id == idx in parent column of edges table
            #    mig_child = np.unique(mig_child)
                
                # Get migration info from nodes list
            #    curr_state = populations[mig_child[0]]
            #    new_state = populations[i]
                
            #    migrant_idx = active_lines.index(mig_child) #change this for ts index
                
                # Update lineage arrays
            #    active_lines[migrant_idx] = idx # name of parent
                
                # Compute event prob
            #    if self.known_ancestral_states:
            #        new_probs = np.zeros(pops)
            #        new_probs[new_state] = 1.0 # event. population
            #        line_state_probs[migrant_idx] = new_probs
            #        event_prob = self.M[curr_state][new_state]
            #    else:
            #        event_prob = 1.0 # pretend as if we don't see migration events
                            
            # Compute prob of no coalescent over time interval
            if not np.isclose(t_elapsed, 0):
                
                if self.known_ancestral_states:
                    
                    # Sum line probs to get total number of lines in each state A
                    A = np.zeros(pops)
                    for probs in line_state_probs: A += probs
                    
                    # Compute prob of no coalescent over time interval
                    pairs = (A * (A-1)) / 2 # number of pairs in each pop
                    lambdas =  pairs * (1/self.Ne) # coal rate in each pop   
                    prob_no_coal = np.exp(-np.sum(lambdas)*t_elapsed)
                
                    # Compute prob of no migration over the time interval
                    #sam = 0
                    #for i in range(pops):
                    #    for z in range(pops):
                    #        sam += (A[i])*(self.M[i][z])
                    #prob_no_mig = np.exp(-sam*t_elapsed)
                    
                    # Compute prob of no reassortment event over the time interval
                    n_segments = np.array(active_segments)
                    # Compute, for every active lineage, the probability that a reassortment event would be observable
                    obs_prob = 1 - (0.5)**(n_segments - 1) 

                    # If pops = 1, this just results in each lineage's obs_prob in a single-column array
                    line_prod = np.array(line_state_probs) * obs_prob[:, np.newaxis]

                    # Summing the reassortment rate across all currently active lineages
                    sum_rate = np.sum(line_prod)  

                    # Scales segment-weighted lineage count by actual per-lineage-per-segment rate parameter
                    prob_no_reassort = np.exp(-sum_rate * reassortment_rate * t_elapsed)

                else: # Unknown ancestral lineage states
                
                    # Integrate lineage prob equations backwards
                    dt_times = list(np.arange(event_time,next_time,self.dt_step)) # integration steps going backwards in time
                    for dt_idx,tx in enumerate(dt_times):
                        
                        # Get time step
                        if (dt_idx+1 < len(dt_times)):
                            dt = dt_times[dt_idx+1] - tx # integration time step
                        else:
                            dt = next_time - tx
    
                        # Should not need to exponentiate transition matrix if dt is small enough
                        expQdt = expm(Q*dt) # exponentiate time-scaled transition rate matrix
    
                        # Update line state probs using Euler integration
                        for ldx,probs in enumerate(line_state_probs):
                            line_state_probs[ldx] = np.matmul(probs,expQdt)
                        
                        # Update total number of lines in each state A
                        A = np.zeros(pops)
                        for probs in line_state_probs: A += probs # sum line probs to get total number of lines in each state
                        
                        # Compute prob of no coalescent over time interval
                        pairs = (A * (A-1)) / 2 # number of pairs in each pop
                        pairs = pairs.clip(min=0) # make sure non are negative
                        lambdas = pairs * (1/self.Ne) # coal rate in each pop
                        prob_no_coal *= np.exp(-np.sum(lambdas)*dt)
                        
                        # Compute prob of no migration over the time interal"
                        #prob_no_mig = 1.0
                        
                        # Compute prob of no reassortment event over the time interval
                        n_segments = np.array(active_segments)
                        obs_prob = 1 - (0.5)**(n_segments - 1)
    
                        line_prod = np.array(line_state_probs) * obs_prob[:, np.newaxis]

                        # Summing the reassortment rate across all currently active lineages
                        sum_rate = np.sum(line_prod)
    
                        prob_no_reassort *= np.exp(-sum_rate * reassortment_rate * dt)

            #log_like += np.log(event_prob) + np.log(prob_no_coal) + np.log(prob_no_mig) + np.log(prob_no_recomb)
            log_like += np.log(event_prob) + np.log(prob_no_coal) + np.log(prob_no_reassort)

        return -log_like


    def _get_line_segment_count(self,line,children,segments):
        
        """
            Compute the number of segments ancestral to the sample and thus eligible to undergo reassortment
        """
        
        line_segments = segments[children == line]
        if len(line_segments) == 0:
            return 0
        all_segs = set()
        for seg_list in line_segments:
            all_segs.update(seg_list)
        return len(all_segs)  
    
    def opt_MLE(self, table_nodes, table_edges):

        """
            Find MLE of single parameter (assumed to be reassortment_rate) using numerical optimization.
            TODO: Generalize to allow for other demographic parameters to be estimataed.
        """

        # Optimize likelihood by minimizing negative log likelihood
        res = minimize_scalar(self.compute_neg_log_like, args=(table_nodes, table_edges), bounds=self.bounds, method='bounded')
        mle = res.x
        
        return mle

if __name__ == '__main__':
       
    #from Espalier.sim import ARGSimulator
    
    # Specify sim params
    samples = 10
    genome_length = 3
    reassortment_rate = 0.01 # NEED TO SET; reassortment rate per lineage per site
    Ne = 100.0  # effective pop sizes
    #M = [[0.0,0.25],[0.25,0.0]]  # migration rate matrix
    M = [[0]]
    
    # Simulation was already run using 20260903_BTB_arg.py
    table_nodes = "arg_w_tables.trees.nodes.csv"
    table_edges = "arg_w_tables.trees.edges.csv"

    # Initialize SCAR model class
    bounds = (0.0,0.1)
    scar_model = SCAR(reassortment_rate,M,Ne,genome_length,bounds=bounds)
    
    # Check numerical optimization for MLE of single param
    mle = scar_model.opt_MLE(table_nodes, table_edges)
    print(mle)
    
    # Check likelihood is valid
    #L = compute_like(ts,**params)
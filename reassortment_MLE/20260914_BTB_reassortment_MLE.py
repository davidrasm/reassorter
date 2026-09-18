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
import sys
import os
import argparse



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
                out_mask = (children == i)          # this sample's one outgoing edge
                out_idx = np.flatnonzero(out_mask)[0]
                active_lines.append(out_idx)
                active_segments.append(len(segments[out_idx]))
                state_probs = np.zeros(pops)
                state_probs[0] = 1.0
                line_state_probs.append(state_probs)

            if 'coalescent' == event_type:
                # incoming: any active lineage whose edge terminates at i
                child_indexes = [k for k, e in enumerate(active_lines) if parents[e] == i]

                coal_probs = np.ones(pops)
                for child_idx in child_indexes:
                    coal_probs *= line_state_probs[child_idx]
                coal_probs = coal_probs / self.Ne
                lambda_sum = sum(coal_probs)
                event_prob = lambda_sum
                parent_probs = coal_probs / lambda_sum

                # outgoing: i's one new edge upward (if any — root has none)
                out_mask = (children == i)
                if np.any(out_mask):
                    out_idx = np.flatnonzero(out_mask)[0]
                    active_lines[child_indexes[0]] = out_idx
                    active_segments[child_indexes[0]] = len(segments[out_idx])
                    line_state_probs[child_indexes[0]] = parent_probs
                    child_indexes.pop(0)
                # Otherwise this is a root: leave child_indexes whole so every # incoming lineage is deleted by the loop below.

                for child_idx in sorted(child_indexes, reverse=True):
                    del active_lines[child_idx]
                    del active_segments[child_idx]
                    del line_state_probs[child_idx]

            if 'recombination' == event_type:
                # incoming: the one active lineage whose edge terminates at i
                child_indexes = [k for k, e in enumerate(active_lines) if parents[e] == i]
                assert len(child_indexes) == 1
                child_idx = child_indexes[0]

                # outgoing: i's edge(s) upward — read parent AND segments straight off them
                out_mask = (children == i)
                out_parents = parents[out_mask]
                out_segments = segments[out_mask]

                if len(out_parents) == 1:
                    # unobservable reassortment — single continuation
                    active_lines[child_idx] = np.flatnonzero(out_mask)[0]
                    active_segments[child_idx] = len(out_segments[0])
                    event_prob = 1.0
                else:
                    assert len(out_parents) == 2, f"Node {i}: expected 1 or 2 parents, got {len(out_parents)}"
                    links = active_segments[child_idx]
                    event_prob = reassortment_rate * (1 - (0.5)**(links - 1))
                    parent_probs = line_state_probs[child_idx]

                    out_idx_0, out_idx_1 = np.flatnonzero(out_mask)

                    active_lines[child_idx] = out_idx_0
                    active_segments[child_idx] = len(out_segments[0])
                    line_state_probs[child_idx] = parent_probs

                    active_lines.append(out_idx_1)
                    active_segments.append(len(out_segments[1]))
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

    # The below function is obselete because the number of segments ancestral to the sample is calculated within the compute_neg_log_likelihood function
    #def _get_line_segment_count(self,line,children,segments):
    #    
    #    """
    #        Compute the number of segments ancestral to the sample and thus eligible to undergo reassortment
    #    """
        
    #    line_segments = segments[children == line]
    #    if len(line_segments) == 0:
    #        return 0
    #    all_segs = set()
    #    for seg_list in line_segments:
    #        all_segs.update(seg_list)
    #    return len(all_segs)  
    
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
    p = argparse.ArgumentParser(
        description="Simulate ancestral recombination graphs (ARGs) under the "
                    "coalescent with recombination or reassortment.")
    p.add_argument("-n", "--num-samples", type=int, default=None,
                    help="number of sampled individuals "
                        "(required via the command line or --config)")
    p.add_argument("-Ne", "--Ne", type=float, default=None,
                    help="effective population size in individuals "
                        "(required via the command line or --config)")
    p.add_argument("--segments", type=int, default=None,
                    help="[reassortment] number of genome segments K")
    p.add_argument("--reassortment-rate", type=float, default=0.0,
                    help="[reassortment] per-lineage reassortment rate")
    p.add_argument("--lower-bound", type=float, required=True,
                    help="lower bound for MLE of reassortment rate")
    p.add_argument("--upper-bound", type=float, required=True,
                    help="upper bound for MLE of reassortment rate")
    p.add_argument("--nodes-csv", required=True,
                    help="path to the .csv with nodes of the ARG")
    p.add_argument("--edges-csv", required=True,
                    help="path to the .csv with edges of the ARG")

    p.set_defaults(simplify=True)

    args = p.parse_args()

    samples = args.num_samples
    genome_length = int(args.segments)
    reassortment_rate = args.reassortment_rate
    Ne = args.Ne  # effective pop sizes
    # Always assuming no migration, for now
    M = [[0]]
    lower_bound = args.lower_bound
    upper_bound = args.upper_bound
    bounds = (lower_bound, upper_bound)
    table_nodes = args.nodes_csv
    table_edges = args.edges_csv
    

    # Initialize SCAR model class
    scar_model = SCAR(reassortment_rate,M,Ne,genome_length,bounds=bounds)
    
    # Check numerical optimization for MLE of single param
    mle = scar_model.opt_MLE(table_nodes, table_edges)
    print(mle)
    
    # Check likelihood is valid
    #L = compute_like(ts,**params)
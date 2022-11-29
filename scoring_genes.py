from network import *
from PF_solutions import *
from utilities import * 
import numpy as np
import time

def gene_substitution(locus_index, itr_gene, solutions, sol_index):
    '''
    Given locus, individual gene, solutions, and solution index,
    perform gene substitution procedure and return chosen genes with the replacement
    '''
    # get the all chosen genes in all the solutions
    sol_chosen_genes = solutions.chosen_genes[sol_index]
    # get the chosen gene at the specified index of the specified solution
    locus_chosen_gene = sol_chosen_genes[locus_index]
    chosen_replaced = list(map(lambda x: x.replace(locus_chosen_gene, itr_gene), sol_chosen_genes))

    return chosen_replaced

def compute_density(chosen_replaced, network):
    '''
    Given chosen genes with replacement, and a network
    density formula: edge counts
    count edges connected among genes in chosen_replaced
    get the number of connections from network (direct neighbors)
    '''
    gene_pairs = [(a, b) for idx, a in enumerate(chosen_replaced) for b in chosen_replaced[idx + 1:]]
    #density = 0
    network_interactions = network.network_interactions
    #subnetwork_interactions = [element for element in gene_pairs if element in network_interactions]
    #subnetwork_interactions = list(set(gene_pairs) & set(network_interactions))
    #subnetwork_interactions = set.intersection(set(gene_pairs), set(network_interactions))
    #subnetwork_interactions = set(gene_pairs).intersection(network_interactions)
    #subnetwork_interactions = set(network_interactions).intersection(gene_pairs)
    density = len(set(gene_pairs).intersection(network_interactions))
    #for gene_pair in gene_pairs:
        #edge_count = network.find_edge(gene_pair[0], gene_pair[1])
        #density = density + edge_count

    return density

def empty_locus_case(locus_index, chosen_replaced, network):
    '''
    Given locus index, chosen genes with replacement, and network 
    Count edges among the chosen genes where the given locus is removed (directed neigbors)
    '''
    # remove the indicated locus completely from the solution
    del chosen_replaced[locus_index]
    gene_pairs = [(a, b) for idx, a in enumerate(chosen_replaced) for b in chosen_replaced[idx + 1:]]
    empty_locus_density = 0
    for gene_pair in gene_pairs:
        edge_count = network.find_edge(gene_pair[0], gene_pair[1])
        empty_locus_density  = empty_locus_density + edge_count

    return empty_locus_density 

def score_gene(locus_index, chosen_replaced, network):
    '''
    Given locus index, chosen genes with replacement, and network,
    compute and return gene score = |edge count(gene sub)-edge count(empty locus case)|
    '''
    density = compute_density(chosen_replaced, network)
    empty_locus_density = empty_locus_case(locus_index, chosen_replaced, network)
    gene_score = np.abs(density - empty_locus_density)

    return gene_score

def get_final_solution(solutions_filename = None, network_filename = None, output_dir=None):
    '''
    Given solutions file name, network file name, and out put directory,
    Compute and return the final scored solution
    Note: use example inputs by default 
    ''' 
    # by default
    if solutions_filename == None:
        loci_candidate_dict,  annotated_candidate_dict = get_loci_candidate_genes()
        # create Prix fixe solutions object
        solutions = PFsolutions(loci_candidate_dict, annotated_candidate_dict)
        # randomly generate chosen genes for each solution (one gene per locus)
        solutions.generate_chosen_genes()
    # solutions given
    else:
         # create Prix fixe solutions object
        chosen_genes, loci_candidate_dict, annotated_candidate_dict = read_in_solutions(chosen_genes_filename = solutions_filename)
        solutions = PFsolutions(loci_candidate_dict, annotated_candidate_dict, chosen_genes)
    
     # create network object
    if network_filename == None:
        network = Network("STRING_network.txt")
    else:
        network = Network(network_filename)

    # get number of solutions = number of all possible combinations of chosen genes
    num_sols = len(solutions.chosen_genes)
    print("Number of input/example solutions: ",num_sols)
    # get the generated chosen genes
    chosen_genes = solutions.chosen_genes
    # iterate through one solution at a time 
    for sol_index in range(1,num_sols):
        # iterate through one locus at a time 
        for locus in solutions.loci_set.keys():
            locus_gene_candidates = solutions.loci_set[locus] 
            locus_chosen_gene = solutions.get_sol_chosen_genes(sol_index)
            locus_itr_genes = list(set(locus_gene_candidates) - set(locus_chosen_gene))
            # calculate score for each candidate gene
            for itr_gene in locus_itr_genes:
                chosen_replaced = gene_substitution(locus, itr_gene, solutions, sol_index)
                gene_score = score_gene(locus, chosen_replaced, network)
                solutions.update_gene_scores(itr_gene, gene_score)
    # average scores across every solution to get final gene scores
    solutions.compute_final_scores()
    # finalize the result dataframe -- add loci annotation
    solutions.finalize_final_sol()
    # output the final solution dataframe
    if output_dir==None:
        output_dir='example_result/final_solution.txt'
    solutions.output_final_sol(output_dir)

    return solutions.final_sol_df
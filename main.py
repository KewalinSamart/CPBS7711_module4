from utilities import *
from PF_solutions import *
from genetic_algorithm import *
from network import *

network = Network("STRING_network.txt")
#chosen_genes_df = pd.read_table('toyexample_solution.txt', delimiter=" ")
#chosen_genes_df.to_csv("toyexample_solution", header=[0,1,2,3,4,5,6,7,8,9,10,11], index=None, sep='\t', mode='a')

chosen_genes, loci_candidate_dict,  annotated_candidate_dict = read_in_solutions('example_solution_cut.txt')
solutions = PFsolutions(loci_candidate_dict,  annotated_candidate_dict, chosen_genes)
#print(solutions.chosen_genes)

loci_indices = sol_locus_tomutate(solutions, percent_mutation=5)
chosen_genes = solutions.chosen_genes
# get candidate genes from the same locus from solutions.loci_set
for locus in range(len(solutions.loci_set.keys())):
    locus_candidates = solutions.loci_set[locus]
    #print(locus_candidates)
    for sol_index in loci_indices[locus]:
        #print(sol_index)
        # pick a gene from the same locus for mutation 
        #print(locus)
        current_chosen_gene = chosen_genes[sol_index-1][locus]
        locus_candiates_removed = list(set(locus_candidates) - set(current_chosen_gene))
        mutated_gene = random.choices(locus_candiates_removed,k=1)[0]
        #print(mutated_gene)
        chosen_genes[sol_index-1][locus] = mutated_gene
# {sol_indices:[[chosen genes at locus0],...,[chosen genes at locus11]]}
# mutated_chosen_genes = chosen_genes
mutated_chosen_genes = dict(zip([k for k in range(1,201)], chosen_genes))
#print(mutated_chosen_genes)
#print(mutated_chosen_genes[1])
#print(mutated_chosen_genes[1][0])
#print(network.network_interactions)


sols_prob_dict = {}
itr = 1
sols_prob = []
sols_density = {} # store solutions' density to compute average density for the generation 
for sol_index in range(1,len(mutated_chosen_genes)+1):
    print(sol_index)
    sol_mutated_chosen = mutated_chosen_genes[sol_index]
    #print(sol_mutated_chosen)
    density = compute_density(sol_mutated_chosen, network)
    sols_density[sol_index] = density
    #print(density)
    #sols_prob.append(density**3)
    sols_prob_dict[sol_index] = density**3
#print(sols_prob_dict)
print("---------sols_prob_dict Done-------------")
sum_cubed_density = sum(sols_prob_dict.values())
sols_prob_dict = {sol_index: cubed_density / sum_cubed_density for sol_index,cubed_density in sols_prob_dict.items()}
#print(sols_prob_dict)
num_solutions = len(sols_prob_dict)
solution_indices = list(sols_prob_dict.keys())
weights_prob = list(sols_prob_dict.values())
# for each iteration, randomly select a pair of sols for num_solutions times
sols_after_mating = {}
for itr_sol in range(1,num_solutions+1): 
    # choose a pair of solutions: list of two sol indices
    random_sol_pair = random.choices(solution_indices,weights=weights_prob,k=2)
    print("random sol pair: ",random_sol_pair)
    # compute average density of the parent subnetworks
    parent_density1 =  compute_density(mutated_chosen_genes[random_sol_pair[0]], network)
    parent_density2 =  compute_density(mutated_chosen_genes[random_sol_pair[1]], network)
    avg_parent_density = (parent_density1+parent_density2)/2
    # sampled_sol_pairs.append(random_sol_pair) # [[sol1,sol10],...,[sol500,sol29]]
    # get mutated genes of all the loci in the two random sols
    mating_res = []
    for locus in range(len(mutated_chosen_genes[1])):
        # choose which solution to get a representative gene from 
        print('locus: ',locus)
        sol_choice = random.choices(random_sol_pair,k=1)
        gene_at_locus = mutated_chosen_genes[int(sol_choice[0])][locus]
        print("chosen gene at the locus from the selected sol: ", gene_at_locus)
        mating_res.append(gene_at_locus)
    sols_after_mating[itr_sol] = mating_res # {sol1:mated_genes,...}
    mating_res_density = compute_density(mating_res, network)
print(sols_after_mating)

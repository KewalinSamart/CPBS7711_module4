from utilities import *
from PF_solutions import *
from genetic_algorithm import *
from network import *

network = Network("STRING_network.txt")
#chosen_genes_df = pd.read_table('example_solutions.txt', delimiter=" ")
#new_df = chosen_genes_df[:5]
chosen_genes, loci_candidate_dict,  annotated_candidate_dict = read_in_solutions('example_solutions.txt',sep=' ')
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
mutated_chosen_genes = dict(zip([k for k in range(1,5001)], chosen_genes))
#print(mutated_chosen_genes)

#sols_prob_dict = {}
itr = 1
sols_prob = []
for sol_index in range(1,len(mutated_chosen_genes)+1):
    print(sol_index)
    sol_mutated_chosen = mutated_chosen_genes[sol_index]
    density = compute_density(sol_mutated_chosen, network)
    print(density)
    sols_prob.append(density)
    #sols_prob_dict[sol_index] = density**3
#sum_cubed_density = sum(sols_prob_dict.values())
#sols_prob_dict = {sol_index: cubed_density / sum_cubed_density for sol_index,cubed_density in sols_prob_dict.items()}

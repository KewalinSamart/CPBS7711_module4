import random 
from PF_solutions import *
from scoring_genes import *
from network import * 
import time

# start with an initial population of 5,000 solutions
# loci_candidate_dict,  annotated_candidate_dict = get_loci_candidate_genes(loci_genes_filename = "toy_loci_set.txt")
# solutions = PFsolutions(loci_candidate_dict, annotated_candidate_dict)
def sol_locus_tomutate(solutions, percent_mutation=5):
    loci_indices = {} 
    num_sols = len(solutions.chosen_genes)
    for locus in range(len(solutions.loci_set.keys())): # locus = 0,1,2,..,12 
        num_to_mutate = int((percent_mutation/100) * num_sols)
        #random.seed(6)
        mut_sol_indices = random.sample([k for k in range(1,num_sols+1)], k=num_to_mutate)
        loci_indices[locus] = mut_sol_indices # value = list of distinct random indices of sols
    return loci_indices

def mutation(solutions, percent_mutation=5):
    # for each locus select 5% of 5,000 solutions to mutate
    # dict = {locus0:[sol_indices -- 5% of 5,000 = 250 distinct indices]}
    # sol_indices are numbers between 1 and 5,000
    loci_indices = sol_locus_tomutate(solutions, percent_mutation)
    #print(loci_indices)
    chosen_genes = solutions.chosen_genes # {1:[chosen_genes],...,N:[chosen_genes]}
    # get candidate genes from the same locus from solutions.loci_set
    for locus in range(len(solutions.loci_set.keys())):
        locus_candidates = solutions.loci_set[locus]
        #print(locus_candidates)
        for sol_index in loci_indices[locus]:
            #print(locus)
            #print(solutions.chosen_genes)
            current_chosen_gene = chosen_genes[sol_index][locus]
            locus_candiates_removed = list(set(locus_candidates) - set(current_chosen_gene))
            mutated_gene = random.choices(locus_candiates_removed,k=1)[0]
            chosen_genes[sol_index][locus] = mutated_gene
    # {sol_indices:[[chosen genes at locus0],...,[chosen genes at locus11]]}
    mutated_chosen_genes = dict(zip([k for k in range(1,len(chosen_genes)+1)], list(chosen_genes.values())))
    return mutated_chosen_genes

def compute_sol_prob(mutated_chosen_genes, network):
    # dict to store normalized cubed density for each subnetwork 
    # where keys: solution numbers, values: normalized probabilities in range(0,1)
    sols_prob_dict = {}
    tic1 = time.perf_counter()
    for sol_index in mutated_chosen_genes.keys():
        sol_mutated_chosen = mutated_chosen_genes[sol_index]
        #print(sol_index)
        tic = time.perf_counter()
        density = compute_density(sol_mutated_chosen, network)
        toc = time.perf_counter()
        print("time spent computing density: ",toc - tic)
        #if density != 0:
        sols_prob_dict[sol_index] = density**3
    sum_cubed_density = sum(sols_prob_dict.values())
    sols_prob_dict = {sol_index: cubed_density / sum_cubed_density for sol_index,cubed_density in sols_prob_dict.items()}
    toc1 = time.perf_counter()
    print("time spent computing probs for the solutions: ",toc1 - tic1)
    return sols_prob_dict

def mating(mutated_chosen_genes, sols_prob_dict, network):
    tic = time.perf_counter()
    solution_indices = list(sols_prob_dict.keys())
    weights_prob = list(sols_prob_dict.values())
    # for each iteration, randomly select a pair of sols for num_solutions times
    sols_after_mating = {}
    mated_density_dict = {}
    for sol in solution_indices: 
        # choose a pair of solutions: list of two sol indices
        random_sol_pair = np.random.choice(solution_indices,size=2,replace=False, p=weights_prob)
        #random_sol_pair = random.sample(solution_indices,weights=weights_prob,k=2)
        # get mutated genes of all the loci in the two random sols
        mating_res = []
        for locus in range(len(mutated_chosen_genes[1])):
            # choose which solution to get a representative gene from 
            choice_index = random.randint(0,1)
            selected_gene = mutated_chosen_genes[random_sol_pair[choice_index]][locus]
            mating_res.append(selected_gene)
        density = compute_density(mating_res, network)
        mated_density_dict[sol] = density
        sols_after_mating[sol] = mating_res # {sol1:mated_genes,...}
        #print("finished mating for solution ",sol)
    toc = time.perf_counter()
    print("time spent mating solutions: ",toc - tic)

    return mated_density_dict, sols_after_mating

def genetic_algorithm(solutions, network, percent_mutation=5):
    prev_gen_avgdensity = 0.0 
    density_list =[prev_gen_avgdensity]
    for i in range(100):
        # the recurring procedures of the mutation and mating steps
        # mutation step
        mutated_chosen_genes = mutation(solutions, percent_mutation)
        #solutions = mutation(solutions, percent_mutation)
        print("Mutation step completed; iteration: ", i)
        # mating step 
        sols_prob_dict = compute_sol_prob(mutated_chosen_genes, network)
        mated_density_dict, sols_after_mating = mating(mutated_chosen_genes, sols_prob_dict, network)
        print("Mating step completed; iteration: ", i)
        print("Parent generation avg density: ",prev_gen_avgdensity)
        # compute average density of the new generation after applying GA
        new_gen_avgdensity = sum(list(mated_density_dict.values()))/len(mated_density_dict.keys())
        density_list.append(new_gen_avgdensity)
        print("New generation avg density: ",new_gen_avgdensity)
        # compare average densities of the new and parent generations
        if ((new_gen_avgdensity - prev_gen_avgdensity) < 0.005*prev_gen_avgdensity):
            break
        else:
            solutions.chosen_genes = sols_after_mating
            prev_gen_avgdensity = new_gen_avgdensity  
            continue
    return density_list, mated_density_dict, solutions
    


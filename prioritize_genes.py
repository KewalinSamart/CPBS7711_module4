from utilities import *
from PF_solutions import *
from genetic_algorithm import *
from network import *
from pval import *
from visualize_finalsol import *
import argparse

def main(output_dir, num_loci=12, lociset_filename='true_disease_loci.txt', network_filename = 'STRING_network.txt', bins = 35, num_lociset=1000, percent_mutation=5, num_solutions=5000, score_cutoff=0.25, num_topscoring=10):
    loci_candidate_dict, annotated_candidate_dict = get_loci_candidate_genes(loci_genes_filename = lociset_filename)
    # initialize true-disease solutions
    solutions = PFsolutions(loci_candidate_dict, annotated_candidate_dict, num_sols = num_solutions)
    # read in network
    network = Network(network_filename)
    avgdensity_list, prev_mated_density_dict, GAoptimized_sols = genetic_algorithm(solutions, network, percent_mutation)
    generate_random_sols(loci_candidate_dict, network, bins = bins, num_lociset=num_lociset, num_sols = num_solutions)
    null_case = generate_null_case(network, percent_mutation, num_lociset, num_solutions)
    pval = calculate_pval_solpop(GAoptimized_sols, null_case)
    # get the final solution with gene scores
    final_sol = get_final_solution(num_loci=num_lociset, network=network, output_dir=output_dir, solutions=GAoptimized_sols)    
    # visualize the result
    finalsol_viz(final_sol, Network(network_filename).network_df, num_loci=num_loci, score_cutoff=score_cutoff)
    # output the final optimized true-disease solutions
    optimized_sols = list(solutions.chosen_genes.values()) # list of lists of chosen genes optimized by GA
    optimized_sols_df = pd.DataFrame(optimized_sols, columns = [i for i in range(len(optimized_sols[0]))])
    optimized_sols_df['density'] = list(prev_mated_density_dict.values())
    optimized_sols_df = optimized_sols_df.sort_values(by = 'density')
    header = [i for i in range(len(optimized_sols[0]))]
    header.append('density')
    optimized_sols_df.to_csv("{}/GAoptimized_sols.txt".format(output_dir), header=header, index=None, sep='\t', mode='a')
    # out put top scoring networks with p-value in the file name
    output_topscoring_networks(GAoptimized_sols, num_topscoring = num_topscoring, pval=pval)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Given a population of solutions for a given set of loci, this script is going to score genes 
    on the loci using the method in Tasan et al. 2015 and visualize the final solution.
    """)
    parser.add_argument("-loci_set", help="(string) solutions file name; set to `toy_loci_set.txt` by default")
    parser.add_argument("-num_loci",help="(int) number of loci; set to `12` by default")
    parser.add_argument("-network", help="txt (tab separated) file containing gene-gene interaction network (undirected; can be weighted/unweighted, but weights will not be used in gene scoring); set to STRING_network.txt by default")
    parser.add_argument("-bins", help="(int) number of bins for gene binning to generate random loci set; set to `35` by default")
    parser.add_argument("-num_lociset", help="(int) number of trials for random set of loci used for null case; set to `1000` by default")
    parser.add_argument("-percent_mutatation",help="(int) percent of mutation of each locus; set to `5` i.e. '5%' by default")
    parser.add_argument("-num_solutions", help=" (int) number of solutions to generate; set to `5000` by default")
    parser.add_argument("-score_cutoff", help="score cutoff for circular-layout (via Networkx) solution visualization; set to `0.25` by default")
    parser.add_argument("top_num", help="number of top scoring networks to output; set to `10` by default")
    parser.add_argument("output_dir", help="(string) path to store final output")


    args = parser.parse_args()
    lociset_filename = args.loci_set
    num_loci = args.num_loci
    network_filename = args.network
    output_dir = args.output_dir
    score_cutoff = args.score_cutoff
    bins = args.bins
    num_lociset = args.num_lociset
    percent_mutation = args.percent_mutation  
    num_solutions = args.num_solutions 
    score_cutoff = args.score_cutoff
    num_topscoring = args.num_topscoring
   

    main(output_dir, num_loci=num_loci,lociset_filename=lociset_filename, network_filename = network_filename, bins = bins, num_lociset=num_lociset, percent_mutation=percent_mutation, num_solutions=num_solutions, score_cutoff=score_cutoff, num_topscoring=num_topscoring)





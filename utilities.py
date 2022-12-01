import pandas as pd 
from network import *

def read_in_solutions(chosen_genes_filename, sep='\t'):
    '''
    Given chosen_genes_filename a file name of solutions, and seperater type
    This function reads in given solutions file in the following format 
    # solutions: 0       1       2      ...
    #            gene3   gene5   gene1  ...
    #            gene4   gene2   gene3  ...
    #            .       .       .      ...
    Returns chosen genes, dictionary of loci-candidate genes {locus:[genes]},
    and annotated dictionary of candidate genes and loci {gene1: locus2,...}
    '''
    chosen_genes_df = pd.read_table(chosen_genes_filename, delimiter=sep)
    # if 'density' in column names then drop the column 
    if 'density' in chosen_genes_df.columns:
        chosen_genes_df.drop('density', inplace=True, axis=1)
    chosen_genes = list(chosen_genes_df.values.tolist())
    loci_candidate_dict = {}
    annotated_candidate_dict = {}
    for locus in chosen_genes_df.columns:
        loci_candidate_dict[int(locus)] = list(set(chosen_genes_df[locus]))
    annotated_candidate_dict = {}
    for key, val in loci_candidate_dict.items():
        for gene in val:
            annotated_candidate_dict[gene] = key

    return chosen_genes, loci_candidate_dict, annotated_candidate_dict

# by default, this program would generate random prix fixe solutions using the gmt file for FA 
def get_loci_candidate_genes(loci_genes_filename = "toy_loci_set.txt"):
    '''
    This function read in txt file with loci and their candidate genes
    input file format:  
    first colunm: loci index e.g. 0, 1, 11
    second column and so on: gene candidate at the specific locus index e.g. AC092291.2
    Returns a dict with keys: loci indices, values: lists of gene candidates
    '''
    myfile = open(loci_genes_filename)
    loci_candidate_dict = {}
    loci_index = 0
    for line in myfile :
        if loci_genes_filename == "toy_loci_set.txt":
            candidate_genes = line.split("\t")[1:]
            candidate_genes = [elem.replace('Locus for ', '') for elem in candidate_genes]
            candidate_genes = [elem.replace('\n', '') for elem in candidate_genes]
        else:
            candidate_genes = line.split("\t")[1:]
        # remove empty strings ("")
        candidate_genes = list(filter(None, candidate_genes))
        loci_candidate_dict[loci_index] =  candidate_genes
        loci_index = loci_index + 1

        annotated_candidate_dict = {}
        for key, val in loci_candidate_dict.items():
            for gene in val:
                annotated_candidate_dict[gene] = key
    return loci_candidate_dict,  annotated_candidate_dict

def set_value(row_num, assigned_value):
    '''
    This is a utility function used to facilitate assigning a distinct color to each locus
    '''
    return assigned_value[row_num]

def output_topscoring_networks(GAoptimized_sols, pval, output_dir, num_topscoring = 10):
    '''
    This function takes a filename of the final optimaized solutions (GA applied) -- txt file
    and num_topscoring (int) a number of topscoring networks to output; 10 by default
    Outputs the top 10 scoring networks (10 highest average densities)
    '''
    topscoring_df =  GAoptimized_sols.tail(10)
    topscoring_df  = topscoring_df.drop('density',axis=1)

    list_nets = list(topscoring_df.values)
    network_df = Network("STRING_network.txt").network_df
    itr = num_topscoring 
    for gene_list in list_nets:
        subnet = network_df.loc[(network_df['gene1'].isin(gene_list)) & (network_df['gene2'].isin(gene_list))]
        subnet.to_csv("{}/topscoring_network{}_{}.txt".format(output_dir,itr,str(pval)), header=['gene1','gene2','weight'], index=None, sep='\t', mode='a')
        itr = itr-1
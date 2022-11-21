import pandas as pd 

def read_in_solutions(chosen_genes_filename, sep='\t'):
    '''
    This function reads in given solutions file in the following format 
    # solutions: 0       1       2      ...
    #            gene3   gene5   gene1  ...
    #            gene4   gene2   gene3  ...
    #            .       .       .      ...
    Returns chosen genes, dictionary of loci-candidate genes {locus:[genes]},
    and annotated dictionary of candidate genes and loci {gene1: locus2,...}
    '''
    chosen_genes_df = pd.read_table(chosen_genes_filename, delimiter=sep)
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
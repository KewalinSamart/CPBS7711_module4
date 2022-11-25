import itertools
import random
import pandas as pd  

class PFsolutions():
    '''
    If no solutions are provided. 
    This class takes information of the default loci set and generates example solutions 
    If a population of solutions is provided, it stores solutions info and gene scores 
    as well as output the final scored solution
    '''
    def __init__(self, loci_candidate_dict, annotated_candidate_dict, chosen_genes = [], num_sols = 5000):
        '''
        Object attributes initialization: loci set, chosen genes, gene scores, 
        final gene scores and final solution dataframe
        '''
        self.loci_set = loci_candidate_dict 
        self.annotated_candidate_dict = annotated_candidate_dict
        if chosen_genes == []:
            chosen_genes = [[] for i in range(num_sols)]
        self.chosen_genes = dict(zip([k for k in range(1,num_sols+1)], chosen_genes))
        #self.chosen_genes = dict(zip([k for k in range(1,len(chosen_genes)+1)], chosen_genes))
        #self.chosen_genes = chosen_genes # empty list by default
        # combine all genes from different loci
        candidate_genes = list(itertools.chain.from_iterable(loci_candidate_dict.values()))
        self.gene_scores = {gene: [0] for gene in candidate_genes} # dict = {gene1:[score_sol1, score_sol2,...],gene2:...}
        self.final_gene_scores = self.gene_scores
        self.final_sol_df = pd.DataFrame()
        
    def generate_chosen_genes(self,iteration=5000):
        '''
        This method generates chosen genes when a population of solutions is not provided
        '''
        chosen_genes = []
        for itr in range(iteration):
            random_genes = []
            for gene_candidates in self.loci_set.values():
                random_gene = random.choice(gene_candidates)
                random_genes.append(random_gene)
            #tp_random_genes = tuple(random_genes)
            chosen_genes.append(random_genes)
        self.chosen_genes = chosen_genes

    def get_sol_chosen_genes(self, sol_index):
        '''
        This method gets the corresponding tuple in chosen_genes list 
        (at position sol_index-1 as sol_index starts with 1)
        '''
        sol_chosen_genes = self.chosen_genes[sol_index-1]

        return list(sol_chosen_genes)
    
    def update_gene_scores(self, gene, score):
        '''
        This method updates list storing gene scores
        '''
        if gene in self.gene_scores.keys():
            self.gene_scores[gene].append(score)
        else:
            self.gene_scores[gene] = [score]
    
    def compute_final_scores(self):
        '''
        This method averages all computed scores for each gene to get the final gene score
        '''
        gene_scores_dict = self.gene_scores
        for key, value in gene_scores_dict.items():
            gene_scores_dict[key] = sum(value)/len(value)
        self.final_gene_scores = gene_scores_dict

    def finalize_final_sol(self):
        '''
        This method finalizes the final output datframe
        '''
        final_scores_df = pd.DataFrame(self.final_gene_scores.items(), columns=['gene','score'])
        loci_candidates_df = pd.DataFrame(self.annotated_candidate_dict.items(), columns=['gene','locus'])
        self.final_sol_df = pd.merge(final_scores_df, loci_candidates_df, on ='gene')
    
    def output_final_sol(self, output_dir='example_result/final_solution.txt'):
        '''
        This method outputs the final scored solution dataframe as a txt file
        '''
        self.final_sol_df.to_csv(output_dir, header=['gene','score','locus'], index=None, sep='\t', mode='a')
        print("The final solution was saved at ", output_dir)
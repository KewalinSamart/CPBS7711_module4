import json
import random 
import numpy as np
import networkx as nx 
import matplotlib.pyplot as plt
from utilities import *
from network import *
import pandas as pd

def generate_loci_colors(num_loci):
    '''
    Given a number of loci, generate a random set of colors (one color for each locus)
    '''
    if num_loci == None: # if not given 
        num_loci = 12
    random_colors =["#"+''.join([random.choice('0123456789ABCDEF') for i in range(6)])for j in range(num_loci)]
    loci_index = np.arange(0,num_loci,1)
    loci_colors = dict(zip(loci_index,random_colors))

    return loci_colors

def generate_networkx_object(final_sol, network_df, loci_colors):
    '''
    Given a final solution (final_sol), network dataframe (network_df),
    and generated loci colors (loci_colors), create a networkx object for visualization and 
    output node attributes
    '''
    genes = list(final_sol.gene.unique())   
    subset_network1 = network_df[network_df.gene1.isin(genes)]
    subset_network2 = subset_network1[subset_network1.gene2.isin(genes)]
    subset_network = subset_network2.dropna()
    G = nx.from_pandas_edgelist(subset_network, 'gene1', 'gene2', "weight")

    node_df = pd.DataFrame(list(G.nodes),columns=['gene'])
    node_attr_df = pd.merge(node_df, final_sol, on ='gene')
    node_attr_df['color'] = node_attr_df['locus'].apply(set_value, args =(loci_colors, ))
    locus_attr = dict(zip(node_attr_df.gene, node_attr_df.locus))
    score_attr = dict(zip(node_attr_df.gene, node_attr_df.score))
    color_attr = dict(zip(node_attr_df.gene, node_attr_df.color))
    nx.set_node_attributes(G, locus_attr , "locus")
    nx.set_node_attributes(G, score_attr , "score")
    nx.set_node_attributes(G, color_attr , "color")

    return G, node_attr_df

def json_cytoscape(G, jsoutput_name="finalsol.json"):
    '''
    Export networkx object to .json for visualization using cytoscape
    ''' 
    json_data = nx.cytoscape_data(G) 
    with open(jsoutput_name, "w") as outfile:
        json.dump(json_data, outfile)

def kamada_kawai_viz(G, node_attr_df):
    '''
    Create and save final solution visualization: kamada_kawai_layout
    '''
    scores = node_attr_df.score
    plt.figure(1, figsize=(150, 80), dpi=40)
    nx.draw(G, node_color=node_attr_df.color, font_size=20, pos=nx.kamada_kawai_layout(G),with_labels = True, node_size=[scores[k]*500 for k in scores],edge_color="grey")
    plt.savefig("finalsol_kkviz.png")   

def circular_viz(G, node_attr_df):
    '''
    Create and save final solution visualization: circular_layout
    '''
    scores = node_attr_df.score
    plt.figure(1, figsize=(50, 25), dpi=50)
    nx.draw(G, node_color=node_attr_df.color, font_size=50, pos=nx.circular_layout(G),with_labels = True, node_size=[scores[k]*1000 for k in scores],edge_color="grey")
    plt.savefig("example_finalsol_ccviz.png")
    
def finalsol_viz(final_sol, network_df, num_loci, score_cutoff):
    '''
    Create and save final solution visualization for both layout types above
    '''
    loci_colors = generate_loci_colors(num_loci)
    G, node_attr_df = generate_networkx_object(final_sol, network_df, loci_colors)
    json_cytoscape(G)
    kamada_kawai_viz(G, node_attr_df)
    if score_cutoff != None:
        final_sol = final_sol[final_sol.score > score_cutoff]
    else:
        final_sol = final_sol[final_sol.score > 0.25]
    G_sub, node_attr_df_sub = generate_networkx_object(final_sol, network_df, loci_colors)
    json_cytoscape(G_sub)
    circular_viz(G_sub, node_attr_df_sub)

def top_scoring_netsviz(solutions_filename, network_filename, num_topscoring = 10, output_dir='topscoring_networks'):
    chosen_genes_df = pd.read_table(solutions_filename, delimiter='\t')
    # get the last 10 rows
    top_scoring_nets = chosen_genes_df.iloc[-num_topscoring:]
    # drop density column
    top_scoring_nets = top_scoring_nets.drop('density', axis=1, inplace=True)
    chosen_genes = list(top_scoring_nets.values.tolist())
    network = Network(network_filename)
    network_interactions = network.network_interactions
    for index in reversed(range(1,len(chosen_genes)+1)):
        gene_set = chosen_genes[index]
        gene_pairs = [(a, b) for idx, a in enumerate(gene_set) for b in gene_set[idx + 1:]]
        topscoring_interactions = set(gene_pairs).intersection(network_interactions)
        # get subset of the network interactions dictionary
        topscoring_interactions = { gene_pair:weight for gene_pair,weight in network_interactions.items() if gene_pair in gene_pairs}
        # turn dict to dataframe
        topscoring_df = pd.DataFrame(topscoring_interactions.items(), columns=['gene1', 'gene2', 'weight'])
        topscoring_df.to_csv("{}/GAoptimized_topscoring_sol{}.txt".format(output_dir,str(index)), header=['gene1', 'gene2', 'weight'], index=None, sep='\t', mode='a')


    

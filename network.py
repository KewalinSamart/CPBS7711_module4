import pandas as pd

class Network():
    '''
    Network class to store and get the given network's properties
    '''
    def __init__(self, network_name, sep = "\t"):
        '''
        read in network and store network dataframe and their gene-gene interactions
        '''
        with open(network_name) as f:
            line = f.readline()
        num_cols = len(line.split())
        if num_cols == 3:
            # network in 3-columns format (weighted)
            network = pd.read_csv(network_name, sep = sep, names = ["gene1","gene2", "weight"],index_col=False)
        elif num_cols == 2:
            # network in 2-columns format (unweighted)
            network = pd.read_csv(network_name, sep = sep, names = ["gene1","gene2"],index_col=False)
            network['weight'] = None 

        # initialize attributes: network data frame, network interactions and degrees
        self.network_df = network
        # set network interactions attribute
        network_list = list(self.network_df.itertuples(index=False, name=None))
        network_interactions = [] 
        for element in network_list:
            strs = list(filter(lambda x : type(x) == str, element))
            network_interactions.append(tuple(sorted(strs)))
        self.network_interactions = network_interactions    
    
    def find_edge(self, gene1, gene2):
        '''
        this method looks into the network interactions and indicate whether there exists an edge 
        connected between the two input genes. if an edge found, return 1, and 0 otherwise
        '''
        key_to_search = tuple(sorted([gene1, gene2]))
        if key_to_search in self.network_interactions:
            edge_found = 1
        else:
            edge_found = 0
        return edge_found
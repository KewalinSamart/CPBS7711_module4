# Disease Gene Prioritization 
### Command
```{r} 
python3 prioritize_genes.py -solutions <solutions file name> -num_loci <number of loci> -network <network file name> -bins <number of bins>, -num_lociset <number of random trials>, -percent_mutatation <chance of mutation> -num_solutions <number of solutions> -score_cutoff <score cutoff> -top_num <number of top scoring networks> output_dir <path to output directory>    
```

### Arguments
- `-solutions` (string) solutions file name; set to `toy_loci_set.txt` by default
- `-num_loci` (int) number of loci; set to `12` by default 
- `-network`  (string) network file name txt (tab separated) file containing gene-gene interaction network (undirected; can be weighted/unweighted, but - weights will not be used in gene scoring); set to `STRING_network.txt` by default
- `-bins` (int) number of bins for gene binning to generate random loci set; set to `35` by default
- `-num_lociset` (int) number of trials for random set of loci used for null case; set to `1000` by default
- `-percent_mutatation` (int) percent of mutation of each locus; set to `5` i.e. 5% by default
- `-num_solutions` (int) number of solutions to generate; set to `5000` by default
- `-score_cutoff` (float) score cutoff for circular-layout (via Networkx) solution visualization; set to `0.25` by default
- `top_num` (int) number of top scoring networks to output; set to `10` by default
- `output_dir` (string) path to store final output 

### Example command
```{r}
python3 prioritize_genes.py true_disease_loci.txt 12 STRING_network.txt 35 1000 5 5000 0.25 10 './path/to/dir'
```

## Details of Inputs/Outputs
### Inputs
1. Undirected un/weighted network file in `.txt` format (`\t` separated) containing gene-gene interactions.
  - `column 1&2`: gene names
  - `column 3`: interaction strengths range from 0 to 1.
  - Example input `STRING_network.txt` –– weighted network:
```{r}
BLOC1S6 BLOC1S3	  0.24
RAB3D   CHML      0.847
MYL7    MYO15A    0.842
...
```
**NOTE:** For unweighted network, network file would be exact the same format as above but without the weight column.

2. 

### Outputs
1. A `.txt` (`\t` separated) file containing genes, their final scores, and the locus that each gene belongs to
- `column 1`: gene names
- `column 2`: final gene score; higher scores indicate better contributions 
- `column 3`: associated locus (locus index)

- Example output `final_solution.txt`
```{r}
gene     score                 locus
PALB2    0.13280434328669868   0
NUPR1    0.02125051082958725   0
SLC5A11  0.016724454415663878  0
CCDC73   0.0                   1
CAPRIN1  0.02042483660130719   1
RCN1     0.05896805896805897   1
...
```
2. Solution visualization in `png` format

Visualization details:
  - Each circle represents a gene
  - Different colors indicate different loci that the genes belong to
  - Circle's size represents their final gene score: a higher score, a bigger circle
  - Edges are gene-gene interaction based on the input network

    2.1 full final solution with no score cutoff applied –– kamada kawai layout
      - Example output `example_result/example_finalsol_kkviz.png`:
      ![alt text](https://github.com/KewalinSamart/CPBS7711_module3/blob/main/example_result/example_finalsol_kkviz.png?raw=true)
    2.2 highest-scores solution with a score cutoff (0.25 by default) –– circular layout
      - Example output `example_result/example_finalsol_ccviz.png`:
      ![alt text](https://github.com/KewalinSamart/CPBS7711_module3/blob/main/example_result/example_finalsol_ccviz.png?raw=true)

3. Final solution subnetwork in `.json` for user to customize the visualization via the interactive tool [`Cytoscape`](https://cytoscape.org/) 
  - Example output `example_finalsol_json.js`:
 ```{r}
 {"data": [], "directed": false, "multigraph": false, "elements": {"nodes": [{"data": {"locus": 8, "score": 0.0358671285918076, "color": "#64678B", 
 "id": "AKAP11", "value": "AKAP11", "name": "AKAP11"}}, {"data": {"locus": 11, "score": 0.1770076779414816, "color": "#D5A612", "id": "SLX4", 
 "value": "SLX4", "name": "SLX4"}}
 ```
 3. Top scoring networks obtimized by the Genetic Algorithm
 

## Installation and Dependencies
- `Python 3.8.3`
- `pandas 1.5.0`
- `numpy 1.23.3`
- `argparse 1.4.0`
- `matplotlib 3.5.1`
- `networkx 2.6`
- `itertools 8.4.0`

 

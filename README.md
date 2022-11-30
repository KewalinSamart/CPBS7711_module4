# Disease Gene Prioritization 
### Command
```{r} 
python3 prioritize_genes.py -solutions <solutions file name> -num_loci <number of loci> -network <network file name> -bins <number of bins>, -num_lociset <number of random trials>, -percent_mutatation <chance of mutation> -num_solutions <number of solutions> -score_cutoff <score cutoff> -top_num <number of top scoring networks> output_dir <path to output directory>    
```

### Arguments
- `-loci_set` (string) solutions file name; set to `toy_loci_set.txt` by default
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

2.gmt `.txt` file (`\t` separated) containing disease-specific loci and their representative disease genes.
```{r}
Fanconi anemia locus 0	Locus for PALB2	NUPR1	CTB-134H23.2	SLC5A11	KIAA0556	CD19	SH2B1	CCDC101	GTF3C1	IL27	ARHGAP17	ERN2	DCTN5	NSMCE1	AQP8	RABEP2	XPO6	ATP2A1	CHP2	BOLA2	KDM8	EIF3C	ATXN2L	LAT	ZKSCAN2	SULT1A1	HS3ST4	EIF3CL	TUFM	NPIPL1	SNX29P2	IL21R	PRKCB	SPNS1	TNRC6A	CACNG3	PLK1	RBBP6	NFATC2IP	APOBR	IL4R	PALB2	SULT1A2	CTD-3203P2.2	GSG1L	SBK1	LCMT1
Fanconi anemia locus 1	Locus for FANCF	CSTF3	FBXO3	SLC17A6	CCDC73	CAPRIN1	RCN1	BDNF	METTL15	CCDC34	EIF3M	LUZP2	BBOX1	CAT	PRRG4	SLC5A12	QSER1	AC103801.2	TCP11L1	SVIP	CD59	NAT10	C11orf91	KCNA4	FIBIN	ARL14EP	ABTB2	LMO2	ELP4	MUC15	DNAJC24	GAS2	LGR4	RP11-17A1.2	MPPED2	KIF18A	FANCF	FSHB	HIPK3	PAX6	DEPDC7	IMMP1L	KIAA1549L	WT1	DCDC5	AC132216.1	ANO5	ANO3	ELF5	EHF	LIN7C
Fanconi anemia locus 2	Locus for RAD51C, BRIP1	CD79B	CACNG1	TANC2	SMG8	RP11-15E18.4	TEX2	YPEL2	RGS9	C17orf72	STRADA	DDX42	TACO1	ICAM2	APOH	PRKCA	FTSJ3	TBX4	DCAF7	GH1	GDPD1	CTD-2510F5.6	METTL2A	MRC2	MAP3K3	PRR11	MED13	C17orf64	TBX2	POLG2	SMURF2	AXIN2	CEP95	INTS2	RAD51C	PPM1E	CA4	CEP112	SMARCD2	C17orf82	USP32	KCNH6	CACNG4	CSH1	RPS6KB1	CTD-2535L24.2	RNFT1	BCAS3	LIMD2	NACA2	RP11-51F16.8	DDX5	APPBP2	SKA2	TRIM37	SCN4A	PTRH2	DHX40	RP11-178C3.1	CCDC47	GNA13	GH2	CSH2	CYB561	HEATR6	VMP1	PSMC5	CSHL1	EFCAB3	TUBD1	CACNG5	BRIP1	PPM1D	AC005544.1	LRRC37A3	CLTC	ERN1	MARCH10	TLK2
```

### Outputs
1. A `.txt` (`\t` separated) file containing genes, their final scores, and the locus that each gene belongs to
- `column 1`: gene names
- `column 2`: final gene score; higher scores indicate better contributions 
- `column 3`: associated locus (locus index)

- Example output `final_solution.txt`
```{r}
gene	  score	 locus
CHP2	  1.5	    0
NSMCE1	0.666666	0
HIPK3	  0.0	    1
DCAF7	  0.75	    2
...
```
2. Solution visualization in `png` format

Visualization details:
  - Each circle represents a gene
  - Different colors indicate different loci that the genes belong to
  - Circle's size represents their final gene score: a higher score, a bigger circle
  - Edges are gene-gene interaction based on the input network

    Highest-scores solution with a score cutoff (0.25 by default) –– circular layout
      - Example output `example_result/example_finalsol_ccviz.png`:
      ![alt text](https://github.com/KewalinSamart/CPBS7711_module3/blob/main/example_result/example_finalsol_ccviz.png?raw=true)

3. Final solution subnetwork in `.json` for user to customize the visualization via the interactive tool [`Cytoscape`](https://cytoscape.org/) 
  - Example output `example_finalsol_json.js`:
 ```{r}
 {"data": [], "directed": false, "multigraph": false, "elements": {"nodes": [{"data": {"locus": 8, "score": 0.0358671285918076, "color": "#64678B", 
 "id": "AKAP11", "value": "AKAP11", "name": "AKAP11"}}, {"data": {"locus": 11, "score": 0.1770076779414816, "color": "#D5A612", "id": "SLX4", 
 "value": "SLX4", "name": "SLX4"}}
 ```
 
 4. Top scoring networks obtimized by the Genetic Algorithm in `.txt` with p-value in the file name
  - Example output `topscoring_network1_0.132.txt`
```{r}
gene1	gene2	weight
UTP3	WDR12	0.856
UTP3	POLR3K	0.158
UTP3	POLR2B	0.268
```
 
## Installation and Dependencies
- `Python 3.8.3`
- `pandas 1.5.0`
- `numpy 1.23.3`
- `argparse 1.4.0`
- `matplotlib 3.5.1`
- `networkx 2.6`
- `itertools 8.4.0`

 

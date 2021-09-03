# Annotation
Scripts to manipulate annotations from Annocript

## Clean contaminants
An R script to identify possible contaminants in the output of Annocript based on the NCBI taxonomy identification.

### Requirements

- Requires the installation of the [Biopython](https://biopython.org/) module in python v3.8;
- Requires installation of [ETE](http://etetoolkit.org/) module (Environment for phylogenetic Tree Exploration) for python v3.8;
- Requires Entrez and NCBITaxa python functions (indicate path in: "load up our functions into memory") ;
- Requires internet access for the Entrez function;

### Usage

'''
Rscript clean_contaminants.R [options]
'''

### Options

'''
   -f CHARACTER, --file=CHARACTER
     Path to the Annocript annotation table
   -o CHARACTER, --out=CHARACTER
     Output prefix
   -m CHARACTER, --mode=CHARACTER
     Search method: locally or online (very slow) [default= locally]
   -s CHARACTER, --fasta=CHARACTER
     Path to Transcripts.fasta sequences file [default= NULL]
   -p CHARACTER, --ORF=CHARACTER
     Path to Annocript_orf_info.fasta file [default= NULL]
   -h, --help
     Show this help message and exit
 '''

### Outputs
- A list with the name of the potentially contaminating transcripts;
- A table with the annotations of the other transcripts in the format provided by Annocript.

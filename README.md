# Annotation
Scripts to manipulate annotations from Annocript

## Clean contaminants
An R script to identify possible contaminants in the output of Annocript based on the NCBI taxonomy identification.

### Requirements

- Requires the installation of the [Biopython](https://biopython.org/) module in python v3.8;
- Requires installation of [ETE](http://etetoolkit.org/) module (Environment for phylogenetic Tree Exploration) for python v3.8;
- Requires Entrez and NCBITaxa python functions (indicate path in: "load up our functions into memory") ;
- Requires internet access for the Entrez function;

### Run it with (last output file name is optional):

### Generates two outputs:
- A list with the name of the potentially contaminating transcripts;
- A table with the annotations of the other transcripts in the format provided by Annocript.

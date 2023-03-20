# Annotation
Scripts to manipulate annotations from Annocript

## Clean Annocript data
An R script to identify possible contaminants in the output of Annocript based on the NCBI taxonomy identification.

### Requirements

- Requires the installation of the [Biopython](https://biopython.org/) module in python v3.8;
- Requires installation of [ETE](http://etetoolkit.org/) module (Environment for phylogenetic Tree Exploration) for python v3.8;
- Requires Entrez and NCBITaxa python functions (indicate path in: "load up our functions into memory") ;
- Requires internet access for the Entrez function;

### Usage

```
Rscript clean_Annocript_data.R [options]
```

### Options

```
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
```

### Outputs
- A list with the name of the potentially contaminating transcripts
- A table of the most frequent contaminants and their frequencies
- A table of transcripts not identified as contaminants
- A .fasta file with the sequences of transcripts not identified as contaminants [optional]
- A .fasta file with the amino acid sequences of the transcript ORFs not identified as contaminants [optional]

## Clean Trinotate data
An R script to identify possible contaminants in the output of Trinotate report.

### Usage

```
Rscript clean_Trinotate_data.R [options]
```

### Options

```
   -t CHARACTER, --trinotate=CHARACTER
     Path to Trinotate report
   -o CHARACTER, --out=CHARACTER
     outprefix name [default= Trinotate]
   -c CHARACTER, --contaminants=CHARACTER
     List of contaminants taxa [default= NULL]
   -a, --assembly
     Logical string (TRUE or FALSE) to define whether to use the assembly data [default= FALSE]
   -f CHARACTER, --fasta=CHARACTER
     Path to Trinity.fasta assembly file [default= Trinity]
   -m CHARACTER, --map=CHARACTER
     Path Trinity.gene_trans_map assembly file [default= Trinity]
   -h, --help
     Show this help message and exit
```

### Outputs
- A table with the name of the potentially contaminating transcripts in the Trinotate format
- A table with the name of non-contaminanting transcripts in the Trinotate format
- A Trinity.fasta file with the sequences of transcripts not identified as contaminants [optional]
- A Trinity.gene_trans_map file of genes not identified as contaminants [optional]

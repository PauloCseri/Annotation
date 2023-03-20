#!/usr/bin/env Rscript --vanila

## ---------------------------
##
## Script name: clean contaminants
##
## Purpose of script: Identify possible contaminants in the output of Annocript 
##  based on the NCBI taxonomy identification. Generates two outputs:
##  - A list with the name of the potentially contaminating transcripts;
##  - A table with the annotations of the other transcripts in the format 
##    provided by Annocript. 
##
## Author: Paulo Cseri Ricardo
##
## Date Created: 2020-06-22
##
## Email: cseri.bio@gmail.com
##
## ---------------------------
##
## Notes: 
##  - Requires the installation of the Biopython module in python v3.8;
##  - Requires installation of ETE module (Environment for phylogenetic Tree 
##    Exploration) n python v3.8;
##  - Requires Entrez and NCBITaxa python functions (indicate path in: "load up 
##  our functions into memory") ;
##  - Requires internet access for the Entrez function;
##  
## ---------------------------
##  
##  Usage: Rscript clean_contaminants.R [options]
##
##
##  Options:
##    -f CHARACTER, --file=CHARACTER
##      Path to the Annocript annotation table
##    -o CHARACTER, --out=CHARACTER
##      Output prefix
##    -m CHARACTER, --mode=CHARACTER
##      Search method: locally or online (very slow) [default= locally]
##    -c CHARACTER,--contaminants=CHARACTER
##      List of contaminants taxa [default= NULL]
##    -s CHARACTER, --fasta=CHARACTER
##      Path to Transcripts.fasta sequences file [default= NULL]
##    -p CHARACTER, --ORF=CHARACTER
##      Path to Annocript_orf_info.fasta file [default= NULL]
##    -h, --help
##      Show this help message and exit
##
##
## ---------------------------

## load up the packages we will need:

cat("\nLoading required packages and functions...\n")

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")

pacman::p_load(optparse,
               reticulate, # reticulate allows use python within R env
               Biostrings,
               stringr) 

## ---------------------------

## load up our (python) functions into env

source_python("/mnt/d/cseri/Documents/Utils_pcr/Entrez.py")
source_python("/mnt/d/cseri/Documents/Utils_pcr/NCBITaxa.py")

## ---------------------------

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Path to annotation file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="Annocript", 
              help="outprefix name [default= %default]", metavar="character"),
  make_option(c("-m", "--mode"), type = "character", default = "locally",
              help = "Search method: locally or online (very slow) [default= %default]",
              metavar = "character"),
  make_option(c("-c", "--contaminants"), type = "character", default = NULL,
              help = "List of contaminants taxa [default= %default]",
              metavar = "character"),
  make_option(c("-s", "--fasta"), type = "character", default = NULL,
              help = "Path to Transcripts.fasta sequences file [default= %default]",
              metavar = "character"),
  make_option(c("-p", "--ORF"), type = "character", default = NULL,
              help = "Path to Annocript_orf_info.fasta file [default= %default]",
              metavar = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (Annocript_filt_ann_out).txt", call.=FALSE)
}


# Define python version (needs v 3.8)
use_python("/usr/bin/python3.8", required = FALSE)    
# Search mode
search = as.character(opt$mode)

## open annotation file
file <- read.table(opt$file, sep = "\t", quote = "", header = TRUE, colClasses = c(rep("character",51)))

## Extract TaxID
tax <- file$Taxonomy %>% str_extract("TaxID.[0-9]*") %>% str_remove("TaxID ")
len <- length(tax)

## New data.frame to receive lines that do not correspond to contaminants
wo_contaminants <- file


## Define COntaminants
if(length(opt$contaminants) == 1){
  Contaminants <- as.character(opt$contaminants)
} else{
  Contaminants <- c("Fungi","Bacteria","Viridiplantae","Acari","Alveolata", "Rhodophyta","Viruses", "Archaea", "metagenomes", "Rhizaria") # taxonomic lineage of contaminants
}

cat("Running:\n")

## Proceed with the analysis of potential contaminants
pb <- txtProgressBar(min = 0, max = len, style = 3) # Create a progress bar

for(i in 1:len){
  taxid <- tax[i]
  if(taxid != "-"){
    
    ## Check which transcripts are in the taxonomic lineage of contaminants (according to NCBI Taxonomy database)
    if(search == "online"){
      lineage <- Entrez(taxid)
    } else if (search == "locally"){
      lineage <- NCBITaxa(strtoi(taxid))
    } else {
      stop("Error: Correctly set the method (locally or online)")
    }
    
    ## Add NAs to contaminants rows
    if(any(lineage %in% Contaminants)){
      wo_contaminants[i,] <- rep(NA,51)
    }  
  }
  setTxtProgressBar(pb,i) # Updates the progress bar
}

wo_contaminants <- wo_contaminants[complete.cases(wo_contaminants),]
contaminants <- file[!file$TranscriptName %in% wo_contaminants$TranscriptName,]

## Count contaminants frequency 
contaminants_TaxID <- contaminants$Taxonomy
ocorrence <- table(contaminants_TaxID)
more.frequent <- ocorrence[ocorrence > 10]

## Save outputs

cat("\nSaving annotations files...\n")
contaminants_transcripts <- contaminants$TranscriptName
no_contaminants_transcripts <- wo_contaminants$TranscriptName
write.table(contaminants_transcripts,paste(opt$out,"_contaminantsID.txt", sep = ""), sep = "\n", col.names = FALSE, row.names = FALSE, quote = FALSE)    # list with contaminants transcripts names
write.table(ocorrence,paste(opt$out,"_contaminants_names.txt", sep = ""),sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(wo_contaminants, paste(opt$out,"_contaminatsfree_annotation.txt", sep = ""), row.names = FALSE, sep = "\t")    # file without contaminants
write.table(more.frequent, paste(opt$out,"_more_frequent_contaminants.txt", sep = ""), row.names = FALSE, sep = "\t")    # file with the frequency of the more frequent contaminants

## save fasta files
if(!is.null(opt$fasta)){
  cat("\nSaving fasta files...\n")
  trans.seq <- readDNAStringSet(filepath = opt$fasta)
  trans.seq.cleanned <- trans.seq[no_contaminants_transcripts]
  writeXStringSet(trans.seq.cleanned,paste(opt$out,"_contaminantsfree.fasta",sep = ""))
}

if(!is.null(opt$ORF)){
  aa.orf <- readAAStringSet(filepath = opt$ORF)
  names(aa.orf) <- gsub(" ].*","",names(aa.orf))
  aa.orf.cleanned <- aa.orf[no_contaminants_transcripts]
  writeXStringSet(aa.orf.cleanned,paste(opt$out,"_orf_contaminantsfree.fasta",sep = ""))
}

cat("\nDone!\n")

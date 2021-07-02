#!/usr/bin/env Rscript

## ---------------------------
##
## Script name: clean_contaminants
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
##  - Requires 3 arguments:
##    @ args[1] Path to the Annocript annotation table
##    @ args[2] Output name
##    @ args[3] Defines the search method, "locally" (by the NCBITaxa function, 
##      recommended) or "online" (by the Entrez function, very slow) 
##
##  Run it with:
##    Rscript clean_contaminants.R Annocript_filt_ann_out.txt outprefix locally
##    or
##    Rscript clean_contaminants.R Annocript_filt_ann_out.txt outprefix online
##
## ---------------------------

## load up the packages we will need:

if(!require(reticulate)){
  install.packages("reticulate")
}
suppressMessages(library(reticulate))  # use python functions within R environment

## ---------------------------

## load up our functions into memory

source_python("/mnt/10T_1/Programas/Utils_pcr/Entrez.py")
source_python("/mnt/10T_1/Programas/Utils_pcr/NCBITaxa.py")

## ---------------------------

args <- commandArgs(trailingOnly = TRUE)
cat("\nStarting program...\n")
use_python("/usr/bin/python3.8", required = FALSE)    # Define python version (needs v 3.8)
search = as.character(args[3])    # Defines where the search will be carried out, locally by the NCBITaxa function (recommended) or online, by the Entrez function (very slow) 

## ###########################
## Open python terminal on R (debug only)
## ###########################
# repl_python()
# exit

## ###########################
## open annotation file
## ###########################
file <- read.table(args[1], sep = "\t", quote = "", header = TRUE, colClasses = c(rep("character",51)))
tax <- file$Taxonomy
len <- length(tax)

## ###########################
## New data.frame to receive lines that do not correspond to contaminants
## ###########################
wo_contaminants <- file[1,]
wo_contaminants <- wo_contaminants[-1,]
contaminants <- file[1,]
contaminants <- contaminants[-1,]

"%ni%" <- Negate("%in%")    # Create %ni% (not in) operator to use in comparisons

## ###########################
## Manually check for contaminants (debug only)
## ###########################
# utax <- sort(unique(tax))
# write.table(utax,"utax.txt")

cat("Running:\n")

## ###########################
# Proceed with the analysis of potential contaminants
## ###########################
pb <- txtProgressBar(min = 0, max = len, style = 3) # Create a progress bar

for(i in 1:len){
  string <- as.character(unlist(tax[i]))
  vector <- unlist(strsplit(string, " "))
  vector.size <- length(vector)
  if(vector.size > 1){
    taxid <- vector[vector.size]
    ## ###########################   
    ## Check which transcripts are in the taxonomic lineage of contaminants (according to NCBI Taxonomy database)
    ## ###########################
    Contaminants <- c("Fungi","Bacteria","Viridiplantae","Acari","Alveolata", "Rhodophyta","Viruses", "Archaea", "metagenomes", "Rhizaria") # taxonomic lineage of contaminants
    delete <- c()
    if(search == "online"){
     lineage <- Entrez(taxid)
    } else {
     lineage <- NCBITaxa(strtoi(taxid))
    }
    
    for(j in 1:length(Contaminants)){
      if(Contaminants[j] %in% lineage){
        delete[j] <- "yes"
      } else {
          delete[j] <- "no"
      }
    }
    ## ###########################
    ## Adds the line file [i,] to the respective data.frame (contaminants or without contaminants) depending on the correspondence with contaminant
    ## ###########################
    if("yes" %ni% delete){
      wo_contaminants[nrow(wo_contaminants) +1,] <- file[i,]
    } else {
     contaminants[nrow(contaminants) +1,] <- file[i,]
    }
  ## ###########################
  ## If there is no taxonomic record (-) add the line file [i,] to the data.frame without contaminants
  ## ###########################
  } else {
      wo_contaminants[nrow(wo_contaminants) +1,] <- file[i,]  
  }
  setTxtProgressBar(pb,i) # Updates the progress bar
}


## ###########################
## Count contaminants frequency 
## ###########################
contaminants_TaxID <- contaminants$Taxonomy
ocorrence <- table(contaminants_TaxID)
more.frequent <- ocorrence[ocorrence > 20]
print(more.frequent)


## ###########################
## Saves the data.frame without contaminants and the list with contaminants transcripts names in files defined by args [2]
## ###########################
cat("\nSaving...")
#cols <- c("TranscriptName","LongOrfLength","ProbToBeNonCoding")
#contaminants$NameFasta <- apply(contaminants[,cols], 1, paste, collapse = "|")
contaminants_transcripts <- contaminants$TranscriptName 
write.table(contaminants_transcripts,paste(args[2],"_contaminantsID.txt", sep = ""), sep = "\n", col.names = FALSE, row.names = FALSE, quote = FALSE)    # list with contaminants transcripts names
write.table(wo_contaminants, paste(args[2],"_contaminantsfree_annotation.txt", sep = ""), row.names = FALSE, sep = "\t")    # file without contaminants
write.table(more.frequent, paste(args[2],"_more_frequent_contaminants.txt", sep = ""), row.names = FALSE, sep = "\t")    # file with the frequency of the more frequent contaminants

cat("\ndone!\n")

#!/usr/bin/env Rscript --vanila

## ---------------------------
##
## Script name: clean contaminants
##
## Purpose of script: Identify possible contaminants in the output of Trinotate 
##  based on the NCBI taxonomy identification. Generates four outputs:
##  - A table with the name of the potentially contaminating transcripts in the
##    Trinotate format;
##  - A table with the name of non-contaminanting transcripts in the Trinotate format.
##  - Two optional outputs from the Trinity assembly data (Trinity.fasta 
##    and Trinity.gene_trans_map) without contaminating transcripts.
##
## Author: Paulo Cseri Ricardo
##
## Date Created: 2023-03-20
##
## Email: cseri.bio@gmail.com
##
## ---------------------------
##  
##  Usage: Rscript clean_Trinotate_data.R [options]
##
##
##  Options:
##    -t CHARACTER, --trinotate=CHARACTER
##      Path to Trinotate report
##    -o CHARACTER, --out=CHARACTER
##      outprefix name [default= Trinotate]
##    -c CHARACTER, --contaminants=CHARACTER
##      List of contaminants taxa [default= NULL]
##    -a, --assembly
##      Logical string (TRUE or FALSE) to define whether to use the assembly data [default= FALSE]
##    -f CHARACTER, --fasta=CHARACTER
##      Path to Trinity.fasta assembly file [default= Trinity]
##    -m CHARACTER, --map=CHARACTER
##      Path Trinity.gene_trans_map assembly file [default= Trinity]
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
               Biostrings,
               stringr,
               dplyr) 

## ---------------------------

option_list = list(
  make_option(c("-t", "--trinotate"), type="character", default=NULL, 
              help="Path to Trinotate report", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="Trinotate", 
              help="outprefix name [default= %default]", metavar="character"),
  make_option(c("-c", "--contaminants"), type = "character", default = NULL,
              help = "List of contaminants taxa [default= %default]",
              metavar = "character"),
  make_option(c("-a","--assembly"), type = "logical", default = FALSE,
              help = "Logical string (TRUE or FALSE) to define whether to use the assembly data [default= %default]"),
  make_option(c("-f", "--fasta"), type = "character", default = NULL,
              help = "Path to Trinity.fasta assembly file [default= %default]",
              metavar = "character"),
  make_option(c("-m", "--map"), type = "character", default = NULL,
              help = "Path Trinity.gene_trans_map assembly file [default= %default]",
              metavar = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$trinotate)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (Trinotate.tsv report)", call.=FALSE)
}


###### Annotation data ###### 
cat("Reading files...\n")

## Open annotation file
# set Trinotate col.names
names <- c("gene_id",
           "transcript_id",
           "sprot_Top_BLASTX_hit",
           "infernal",
           "prot_id",
           "prot_coords",
           "sprot_Top_BLASTP_hit",
           "Pfam",
           "SignalP",
           "TmHMM",
           "eggnog",
           "Kegg",
           "gene_ontology_BLASTX",
           "gene_ontology_BLASTP",
           "gene_ontology_Pfam",
           "transcript",
           "peptide")

# Open Trinotate file
Trinotate <- read.table(opt$trinotate, 
                        skip = 3, 
                        header = FALSE, 
                        sep = "\t",
                        fill = TRUE,
                        quote = "",
                        col.names = names)


## set Contaminants
if(length(opt$contaminants) == 1){
  contaminants <- as.character(opt$contaminants)
} else{
  contaminants <- c("Fungi",
                  "Bacteria",
                  "Viridiplantae",
                  "Acari",
                  "Alveolata",
                  "Rhodophyta",
                  "Viruses",
                  "Archaea", 
                  "metagenomes", 
                  "Rhizaria") # taxonomic lineage of contaminants
}

## create vector of genes that showed blast hits with contaminants
cat("Running:\n")

n.Trinotate <- nrow(Trinotate) 
contaminant.genes <- NULL # contaminant genes vector
contaminant <- NULL

pb <- txtProgressBar(min = 0, max = n.Trinotate, style = 3) # Create a progress bar
for(i in 1:n.Trinotate){
  char <- Trinotate$sprot_Top_BLASTX_hit[i] %>% str_split(";\\^") %>% unlist() %>%
    nth(2) %>% str_split("; ") %>% unlist() 
  test <- char %in% contaminants %>% any()
  if(test){
    contaminant.genes <- append(contaminant.genes,Trinotate$gene_id[i])
    contaminant <- append(contaminant,char[char %in% contaminants])
  }
  setTxtProgressBar(pb,i) # Updates the progress bar
}

# contingency table
cont <- data.frame(gene_id=contaminant.genes,
                   contaminant=contaminant)
cat("\n \nContaminating sequences:\n")
print(table(cont$contaminant))

# filter out contaminant gene annotations
contaminant.genes <- unique(contaminant.genes)
contaminants.annot <- Trinotate %>% filter(gene_id %in% contaminant.genes)

# filter out non-contaminant gene annotations
contaminants.free <- Trinotate %>% filter(!gene_id %in% contaminant.genes)


## Save annotation outputs
cat("\nSaving annotations files...\n")

header <- "Use of uninitialized value $blast_type in concatenation (.) or string at /mnt/10T_1/Programas/Trinotate-Trinotate-v4.0.0/Trinotate line 200.
-REPORT being generated.
CMD: /mnt/10T_1/Programas/Trinotate-Trinotate-v4.0.0/util/Trinotate_report_writer.pl --sqlite /mnt/10T_2/larissa/annotation/Trinotate_03_2023.sqlite -E 1e-5 --pfam_cutoff DNC
#gene_id        transcript_id   sprot_Top_BLASTX_hit    infernal        prot_id prot_coords     sprot_Top_BLASTP_hit    Pfam    SignalP TmHMM   eggnog  Kegg    gene_ontology_BLASTX    gene_ontology_BLASTP    gene_ontology_Pfam      transcript      peptide"

# contamination free annotation
write(header,paste("contaminants_free_annot_", opt$out, "_", Sys.Date(),".tsv",sep = ""))
write.table(contaminants.free, 
            paste("contaminants_free_annot_", opt$out, "_", Sys.Date(),".tsv",sep = ""),
            append = TRUE,
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

# Annotation from contaminant transcripts
write(header,paste("contaminants_annot_", opt$out, "_", Sys.Date(),".tsv", sep = ""))
write.table(contaminants.annot,
            paste("contaminants_annot_", opt$out, "_", Sys.Date(),".tsv",sep = ""),
            append = TRUE,
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

######  Assembly data ######
if(opt$assembly){
  ## Open assembly files
  # open map file
  maps <- read.table(opt$map,
                     sep = "\t", col.names = c("gene_id","transcript_id"))
  # open fata file
  fasta <- readDNAStringSet(filepath = opt$fasta)
  
  ## select non-contaminant data from assembly data
  contaminants.free.maps <- maps %>% filter(! gene_id %in% contaminant.genes)
  contaminant.free.transcripts <- contaminants.free.maps %>% pull(transcript_id)
  names.fasta <- names(fasta)
  pos <- gsub(" len=.*", "",names.fasta) %in% contaminant.free.transcripts
  contaminants.free.fasta <- fasta[pos]
  
  ## Saving assembly files  
  cat("\nSaving assembly files...\n")
  
  # save map file
  write.table(contaminants.free.maps,
    paste("contaminants_free_",
          opt$out,
          "_",
          Sys.Date(),
          ".gene_trans_map",sep = ""),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t")
  
  # save fasta file
  writeXStringSet(contaminants.free.fasta,
                  paste("contaminants_free_",
                        opt$out, 
                        "_",
                        Sys.Date(),
                        ".fasta",sep = ""))
}
cat("\nDone!\n")

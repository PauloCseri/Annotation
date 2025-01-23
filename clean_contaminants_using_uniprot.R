#!/usr/bin/env Rscript

## ---------------------------
## Script name: clean_contaminants_using_uniprot.R
##
## Purpose: Identify possible contaminants in the output of Blastx and Blastp against
## UniRef90 based on NCBI taxonomy identification. Generates clean
## outputs for further analysis.
##
## Author: Paulo Cseri Ricardo
## Email: cseri.bio@gmail.com
##
## Last Updated: 2025-01-17
## ---------------------------

# Load required packages
cat("\nLoading required packages and functions...\n")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(optparse, Biostrings, stringr, dplyr, tidyr, writexl, data.table, progress)

# Define options for script usage
option_list <- list(
  make_option(c("-b", "--blastxout"), type = "character", default = NULL, 
              help = "Path to Blastx out file", metavar = "character"),
  make_option(c("-p", "--blastpout"), type = "character", default = NULL, 
              help = "Path to Blastp out file", metavar = "character"),
  make_option(c("-t", "--trinotate"), type = "character", default = NULL, 
              help = "Path to Trinotate out file", metavar = "character"),
  make_option(c("-u", "--uniprotx"), type = "character", default = NULL, 
              help = "Path to Blastx Uniprot ID mapping table", metavar = "character"),
  make_option(c("-v", "--uniprotp"), type = "character", default = NULL, 
              help = "Path to Blastp Uniprot ID mapping table", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "Trinotate", 
              help = "Output prefix name [default: %default]", metavar = "character"),
  make_option(c("-c", "--contaminants"), type = "character", default = NULL, 
              help = "List of contaminants taxa [default: NULL]", metavar = "character"),
  make_option(c("-a", "--assembly"), type = "logical", default = FALSE, 
              help = "Logical flag to define use of assembly data [default: FALSE]"),
  make_option(c("-f", "--fasta"), type = "character", default = NULL, 
              help = "Path to Trinity.fasta assembly file", metavar = "character"),
  make_option(c("-m", "--map"), type = "character", default = NULL, 
              help = "Path to Trinity.gene_trans_map assembly file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
required_args <- c("blastxout", "blastpout", "trinotate", "uniprotx", "uniprotp")
missing_args <- required_args[!sapply(required_args, function(arg) !is.null(opt[[arg]]))]
if (length(missing_args) > 0) {
  print_help(opt_parser)
  stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")), call. = FALSE)
}

# Load input files
cat("\nReading input files...\n")
blastx_out <- fread(opt$blastxout, header = FALSE, sep = "\t", fill = TRUE)
blastp_out <- fread(opt$blastpout, header = FALSE, sep = "\t", fill = TRUE)
uniprotx <- fread(opt$uniprotx, header = TRUE, sep = "\t", fill = TRUE)
uniprotp <- fread(opt$uniprotp, header = TRUE, sep = "\t", fill = TRUE)

# Load Trinotate file
trinotate_columns <- c(
  "gene_id", "transcript_id", "sprot_Top_BLASTX_hit", "infernal", "prot_id",
  "prot_coords", "sprot_Top_BLASTP_hit", "Pfam", "SignalP", "TmHMM", "eggnog",
  "Kegg", "gene_ontology_BLASTX", "gene_ontology_BLASTP", "gene_ontology_Pfam",
  "transcript", "peptide"
)
trinotate <- fread(opt$trinotate, skip = 4, sep = "\t", fill = TRUE, col.names = trinotate_columns)

# Clean lineage strings
clean_lineage <- function(lineage_column) {
  lineage <- as.character(lineage_column) 
  lineage <- gsub(" \\(.*?\\)", "", lineage) # Remove parenthetical content
  lineage <- gsub(",", ";", lineage) # Replace commas with semicolons
  return(lineage)
}

uniprotx$Taxonomic.lineage <- clean_lineage(uniprotx$`Taxonomic lineage`)
uniprotp$Taxonomic.lineage <- clean_lineage(uniprotp$`Taxonomic lineage`)

# Merge data
merge_with_uniprot <- function(blast_out, uniprot_table) {
  uniprot_table <- uniprot_table[!duplicated(uniprot_table$From),]
  merged <- base::merge(blast_out, uniprot_table, by.x = "V2", by.y = "From", all.x = TRUE)
  merged[is.na(merged)] <- "."
  merged[merged == ""] <- "."
  return(merged)
}
merged_blastx <- merge_with_uniprot(blastx_out, uniprotx)
merged_blastp <- merge_with_uniprot(blastp_out, uniprotp)

# Filter contaminants
cat("\nFiltering contaminants...\n")
if (is.null(opt$contaminants)) {
  contaminants <- c("Fungi", "Bacteria", "Viridiplantae", "Acari", "Alveolata",
                    "Rhodophyta", "Viruses", "Archaea", "metagenomes", "Rhizaria")
} else {
  contaminants <- unlist(strsplit(opt$contaminants, ","))
}

filter_contaminants <- function(merged_data, contaminants) {
  progress_bar <- progress::progress_bar$new(
    total = nrow(merged_data),
    format = "[:bar] :percent Remaining: :eta"
  )
  filtered <- merged_data %>%
    mutate(is_contaminant = sapply(Taxonomic.lineage, function(lineage) {
      progress_bar$tick()
      any(unlist(strsplit(lineage, "; ")) %in% contaminants)
    })) %>%
    filter(!is_contaminant) %>%
    select(-is_contaminant)
  return(filtered)
}
cat("\nFiltering UniRef blastx data...\n")
contaminant_free_blastx <- filter_contaminants(merged_blastx, contaminants)
cat("\nFiltering UniRef blastp data...\n")
contaminant_free_blastp <- filter_contaminants(merged_blastp, contaminants)

# Convert GO annotations to trinotate format
cat("\nMerging and editing data...\n")
convert_go_annotations <- function(column_vector, ontology) {
  # Apply the conversion to each element in the vector
  result <- sapply(column_vector, function(x) {
    if (x == ".") {
      return(x) # Keep the "." values unchanged
    } else {
      # Remove line breaks (\n) and leading/trailing spaces
      x <- gsub("\n", " ", x)           # Remove line breaks
      x <- trimws(x, which = "both")    # Remove leading and trailing spaces
      # Split the original string by ";"
      annotations <- strsplit(x, ";\\s*")[[1]]
      # Convert each annotation to the new format
      converted_annotations <- sapply(annotations, function(annotation) {
        # Extract the GO ID and description
        match <- regmatches(annotation, regexec("^(.*) \\[(GO:\\d+)\\]$", annotation))
        if (length(match[[1]]) > 0) {
          description <- match[[1]][2] # The term (e.g., "DNA topoisomerase binding")
          go_id <- match[[1]][3]       # The GO ID (e.g., "GO:0044547")
          return(paste(go_id, ontology, description, sep = "^"))
        } else {
          # Return the original annotation if the pattern is not recognized
          return(annotation)
        }
      })
      # Join the converted annotations with "`"
      return(paste(converted_annotations, collapse = "`"))
    }
  })
  
  return(result)
}

# Convert in blastx data
contaminant_free_blastx$`Gene Ontology (biological process)` <- convert_go_annotations(contaminant_free_blastx$`Gene Ontology (biological process)`, "biological_process") 
contaminant_free_blastx$`Gene Ontology (cellular component)` <- convert_go_annotations(contaminant_free_blastx$`Gene Ontology (cellular component)`, "cellular_component")
contaminant_free_blastx$`Gene Ontology (molecular function)` <- convert_go_annotations(contaminant_free_blastx$`Gene Ontology (molecular function)`, "molecular_function")

# Convert in blastp data
contaminant_free_blastp$`Gene Ontology (biological process)` <- convert_go_annotations(contaminant_free_blastp$`Gene Ontology (biological process)`, "biological_process") 
contaminant_free_blastp$`Gene Ontology (cellular component)` <- convert_go_annotations(contaminant_free_blastp$`Gene Ontology (cellular component)`, "cellular_component")
contaminant_free_blastp$`Gene Ontology (molecular function)` <- convert_go_annotations(contaminant_free_blastp$`Gene Ontology (molecular function)`, "molecular_function")


# Unite and correct wrong patterns
clean_strings <- function(column_vector) {
  # Replace ".`.`." with "."
  column_vector <- gsub("\\.\\`\\.\\.", ".", column_vector)
  # Replace ".`" with ""
  column_vector <- gsub("\\.\\`", "", column_vector)
  # Replace "`." with ""
  column_vector <- gsub("\\`\\.", "", column_vector)
  return(column_vector)
}

contaminant_free_blastx <- contaminant_free_blastx %>% unite(
  gene_ontology_BLASTX_Uniprot, `Gene Ontology (cellular component)`, 
  `Gene Ontology (molecular function)`,`Gene Ontology (biological process)`, sep = "`",
  remove = FALSE) %>% mutate(gene_ontology_BLASTX_Uniprot = 
                               clean_strings(gene_ontology_BLASTX_Uniprot))

contaminant_free_blastp <- contaminant_free_blastp %>% unite(
  gene_ontology_BLASTP_Uniprot, `Gene Ontology (cellular component)`, 
  `Gene Ontology (molecular function)`,`Gene Ontology (biological process)`, sep = "`",
  remove = FALSE) %>% mutate(gene_ontology_BLASTP_Uniprot = 
                               clean_strings(gene_ontology_BLASTP_Uniprot))


# Combine columns to match to Trtinotate pattern
combine_columns <- function(df, v2, entry_name, v7, v8, v9, v10, v3, v11, protein_names, taxonomic_lineage, new_column_name) {
  df %>%
    mutate(
      !!new_column_name := paste(
        .[[v2]],
        .[[entry_name]],
        paste0("Q:", .[[v7]], "-", .[[v8]]),
        paste0("H:", .[[v9]], "-", .[[v10]]),
        paste0(.[[v3]], "%ID"),
        paste0("E:", .[[v11]]),
        .[[protein_names]],
        .[[taxonomic_lineage]],
        sep = "^"
      )
    )
}


contaminant_free_blastx <- combine_columns(contaminant_free_blastx, "V2",
                                           "Entry Name", "V7","V8", "V9", "V10",
                                           "V3", "V11", "Protein names", 
                                           "Taxonomic.lineage", "uniprot_Top_BLASTX_hit")
contaminant_free_blastx <- contaminant_free_blastx %>% rename(transcript_id = V1) %>% 
  select(uniprot_Top_BLASTX_hit, transcript_id, gene_ontology_BLASTX_Uniprot)



contaminant_free_blastp <- combine_columns(contaminant_free_blastp, "V2",
                                           "Entry Name", "V7","V8", "V9", "V10",
                                           "V3", "V11", "Protein names", 
                                           "Taxonomic.lineage", "uniprot_Top_BLASTP_hit")
contaminant_free_blastp <- contaminant_free_blastp %>% rename(prot_id = V1) %>%
  select(uniprot_Top_BLASTP_hit, prot_id, gene_ontology_BLASTP_Uniprot)


# Merge annotations
merge_annotations <- function(trinotate_df, uniprot, blast, molecule) {
  # Perform a join based on the transcript_id column
  merged <- dplyr::left_join(trinotate_df, uniprot, by = paste0(molecule,"_id"))
  
  # Keep only the desired observations
  keep_rows <- merged %>%
    filter(!is.na(!!sym(paste0("gene_ontology_", blast,"_Uniprot"))) | 
             (is.na(!!sym(paste0("gene_ontology_", blast,"_Uniprot"))) & sprot_Top_BLASTX_hit == "." & sprot_Top_BLASTP_hit == "."))
  
  keep_rows[is.na(keep_rows)] <- "."
  
  return(keep_rows)
}



annotations <- merge_annotations(trinotate, contaminant_free_blastx, "BLASTX", "transcript")
annotations <- merge_annotations(annotations, contaminant_free_blastp, "BLASTP", "prot")
annotations <- annotations %>% select(gene_id, transcript_id, sprot_Top_BLASTX_hit,
                                      uniprot_Top_BLASTX_hit, infernal, prot_id, prot_coords,
                                      sprot_Top_BLASTP_hit, uniprot_Top_BLASTP_hit, Pfam,
                                      SignalP, TmHMM, eggnog, Kegg, gene_ontology_BLASTX, gene_ontology_BLASTP,
                                      gene_ontology_Pfam, gene_ontology_BLASTX_Uniprot, gene_ontology_BLASTP_Uniprot,
                                      transcript, peptide)

# Filter swiss-prot contaminants
filter_contaminants2 <- function(uniref_filtered, contaminants, blast) {
  progress_bar <- progress::progress_bar$new(
    total = nrow(uniref_filtered),
    format = "[:bar] :percent Remaining: :eta"
  )
  
  filtered <- uniref_filtered %>%
    mutate(is_contaminant = mapply(function(uniprot_hit, sprot_hit) {
      progress_bar$tick()
      
      # If uniprot or swiss-prot contains only ".", keep the row
      if (sprot_hit == ".") {
        return(FALSE)
      }
      
      # Check if the last element of the uniprot string is "."
      uniprot_unknow <- FALSE
      uniprot_lineage <- sub(".*\\^", "", uniprot_hit)
      uniprot_parts <- unlist(strsplit(uniprot_lineage, "; "))
      if (tail(uniprot_parts, 1) == ".") {
        uniprot_unknow <- TRUE
      } else {
        return(FALSE)
      }
      
      # Check if swiss-prot contains contaminants
      if(uniprot_unknow){
        sprot_lineage <- sub(".*\\^", "", sprot_hit)
        sprot_parts <- unlist(strsplit(sprot_lineage, "; "))
        if (any(sprot_parts %in% contaminants)) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      }
    }, uniref_filtered[[sym(paste0("uniprot_Top_",blast,"_hit"))]], 
    uniref_filtered[[sym(paste0("sprot_Top_",blast,"_hit"))]])) %>%
    filter(!is_contaminant) %>%
    select(-is_contaminant)
  
  return(filtered)
}

cat("\nFiltering UniRef blastx data...\n")
annotations <- filter_contaminants2(annotations, contaminants, "BLASTX")
cat("\nFiltering UniRef blastp data...\n")
annotations <- filter_contaminants2(annotations, contaminants, "BLASTP")


# get contaminants annotations
non_contaminants <- annotations$gene_id
contaminant.genes <- trinotate$gene_id[!trinotate$gene_id %in% non_contaminants]
contaminants_annotations <- trinotate %>% filter(gene_id %in% contaminant.genes)

# Output results
cat("\nWriting output files...\n")
fwrite(annotations, paste0("contaminants_free_annot_", opt$out, "_", Sys.Date(),".xls"), sep = "\t")
fwrite(annotations, paste0("contaminants_free_annot_", opt$out, "_", Sys.Date(),".tsv"), sep = "\t")
fwrite(contaminants_annotations, paste0("contaminants_annot_", opt$out, "_", Sys.Date(),".tsv"), sep = "\t")


######  Assembly data ######
if(opt$assembly){
  ## Open assembly files
  # open map file
  maps <- read.table(opt$map,
                     sep = "\t", col.names = c("gene_id","transcript_id"))
  # open fata file
  fasta <- readDNAStringSet(filepath = opt$fasta)
  
  ## select non-contaminant data from assembly data
  contaminants.free.maps <- maps %>% filter(gene_id %in% non_contaminants)
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
cat("\nProcessing completed successfully!\n")


####################################################################################################
# SCRIPT : RNA-seq SLA - Extraction, annotation & normalisation des gènes d’intérêt (niveau gène)
# Auteur : Mickael Coquerelle
# Date : 09/07/2025
# Objectif : Générer une matrice propre et annotée des 56 gènes SLA à partir d’un comptage génique.
####################################################################################################

# Chargement des bibliothèques
library(data.table)
library(biomaRt)
library(tidyverse)


# Connexion à Ensembl :
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extraction des correspondances Gene ID / Gene Name
mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),mart = ensembl, useCache = FALSE)

# Renommage des colonnes
colnames(mapping) <- c("Geneid", "gene_name")
write.table(mapping, file = "ENSG_to_Name2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Fichiers d'entrée
count_file <- "all-counts.tsv"
mapping_file <- "ENSG_to_Name.tsv"

# Liste des 56 gènes SLA d’intérêt du panel ciblé
genes_list <- c(
  "ALS2", "ANG", "ANO6", "APOE", "APP", "C21orf2", "C9orf72", "CAPRIN1", "CCNF", "CHCHD10",
  "CHMP2B", "CSF1R", "DAO", "DCTN1", "DPYSL3", "ELP3", "EPHA4", "ERBB4", "EWSR1", "FAS",
  "FIG4", "FUS", "GC", "GRN", "HNRNPA1", "HNRNPA2B1", "KIF5A", "MATR3", "MOBP", "NEFH",
  "NEK1", "OPTN", "PFN1", "PRPH", "SARM1", "SOD1", "SPG11", "SQSTM1", "TARDBP", "TIA1",
  "TREM2", "TRPM7", "TUBA4A", "UBQLN2", "UNC13A", "USP14", "USP7", "VAPB", "VCP", "SETX",
  "SIGMAR1", "ANGEL2", "TBK1", "ZFHX3", "ZNF512"
)
all_counts <- fread(count_file)

# Nettoyage des noms de colonnes BAM
bam_cols <- colnames(all_counts)[7:ncol(all_counts)]
clean_names <- sub("\\.bam$", "", basename(bam_cols))
colnames(all_counts)[7:ncol(all_counts)] <- clean_names # remplacement par les noms de patients nettoyés

# Annotation avec les noms de gènes
mapping <- fread(mapping_file, col.names = c("Geneid", "gene_name"))
all_counts_annotated <- left_join(all_counts, mapping, by = "Geneid")

all_counts_sla <- all_counts_annotated %>% filter(tolower(gene_name) %in% tolower(genes_list)) %>% select(Geneid,gene_name,everything()) %>%
  select(-c(Chr,Start,End,Strand))

col_info <- tibble(full_name = colnames(all_counts)[7:ncol(all_counts)])
#view(col_info)

genes_found <- unique(tolower(all_counts_sla$gene_name))
#View(genes_found)
genes_missing <- setdiff(tolower(genes_list), genes_found)
if (length(genes_missing) > 0) {
  warning("Gènes SLA non retrouvés : ", paste(genes_missing, collapse = ", "))
}

# Sauvegarde du tableau filtré et annoté
write.table(all_counts_sla, file = "genes_SLA_56_filtered.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

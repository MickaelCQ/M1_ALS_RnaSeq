############################################################################################################
# SCRIPT : RNA-seq SLA - Prétraitement, annotation & extraction des gènes d’intérêt
# Auteur : Mickael Coquerelle
# Date : 28/05/2025
# Description :La première partie du script (Etape 1 à 6) 
#              Ce script traite une matrice de comptage RNA-seq issue de featureCounts
#              j'ajoute les noms de gènes  à partir d’un fichier GTF GRCh37 des correspondances avec les "ENSG", 
#              puis une extraction des 56 gènes d’intérêt SLA. avec une mise en forme, propre. du tsv de sortie,
#              Pour manipulations futures.

#              La deuxieme partie du script (Etape 7 et plus), consiste à formaliser les données en data.frame 
#              cohérent pour effectuer les différents scénario de normalisation 

#              La troisième partie du script (Etape N ...), 
############################################################################################################


############################################################################################################
#                                 Chargement des librairies nécessaires 
############################################################################################################
library(tidyverse) # Wickham H. (2017) "Tidy Data". Journal of Statistical Software.
library(edgeR)    # Robinson MD et al. (2010), Bioinformatics. Analyse des données RNA-seq (comptage, normalisation)
library(biomaRt)  # Durinck S. et al. (2009), Nature Protocols. Accès à Ensembl depuis R pour avoir les annotations.
library(DESeq2)   # Love MI et al. (2014), Genome Biology. 
library(ggplot2)
library(data.table)

############################################################################################################
#                 Étape 1 : Extraction des noms de gènes depuis un fichier GTF 
############################################################################################################
# Fichier GTF d’annotation Ensembl GRCh37
gtf_file <- "~/Téléchargements/Homo_sapiens.GRCh37.87.chr.gtf.gz"

# On extrait uniquement les lignes contenant des annotations de type "gene"
gtf <- fread(cmd = paste("zgrep -P '\tgene\t' ", gtf_file), header = FALSE, sep = "\t", data.table = FALSE)

# On récupère les attributs pour extraire gene_id et gene_name
attr_field <- gtf[, 9]

extract_attr <- function(attr_str, attr_name) {
  pattern <- paste0(attr_name, ' "([^"]+)"')
  m <- regmatches(attr_str, regexec(pattern, attr_str))
  sapply(m, function(x) if(length(x) > 1) x[2] else NA)
}

gene_id   <- extract_attr(attr_field, "gene_id")
gene_name <- extract_attr(attr_field, "gene_name")

# Table de correspondance ENSG match avec les noms de gènes
mapping <- data.frame(gene_id = gene_id, gene_name = gene_name, stringsAsFactors = FALSE)

# Sauvegarde locale
write.table(mapping, file = "ENSG_to_Name.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

############################################################################################
#                Étape 2 : Chargement de la matrice de comptage featureCounts              #
############################################################################################
df <- read.delim("all-counts.tsv", sep = '\t', header = TRUE, skip = 1, check.names = FALSE)

# Extraction des colonnes d’expression (issues des fichiers BAM)
bam_cols <- grep("^output/STAR/bam/", colnames(df), value = TRUE)
counts   <- df[, bam_cols]

# Extraction des informations sur les gènes
gene_info <- df[, c("Geneid", "Chr", "Start", "End", "Strand", "Length")]

# Conversion en matrice numérique
counts_mat <- as.matrix(counts)
mode(counts_mat) <- "numeric"

# Longueurs des gènes
lengths <- gene_info$Length
df$Chr  <- sapply(strsplit(as.character(df$Chr), ";"), `[`, 1)

# Nettoyage des noms de colonnes BAM pour que ce soit plus lisible pr les ggplot 
clean_colnames <- function(x) {
  x <- basename(x)                # garde juste le nom du fichier
  x <- gsub("\\.bam$", "", x)     # retire l'extension
  return(x)
}
new_bam_names <- clean_colnames(bam_cols)
colnames(df)[match(bam_cols, colnames(df))] <- new_bam_names
############################################################################################################
#                       Étape 3 : Filtrage des gènes exprimés au moins une fois                            #
############################################################################################################

expressed_genes <- rowSums(counts_mat) > 1
counts_mat_filt <- counts_mat[expressed_genes, ]
lengths_filt    <- lengths[expressed_genes]

############################################################################################################
#                       Étape 4 : Annotation des gènes avec les noms officiels (HGNC)                      #
############################################################################################################
Correspondance <- read.delim("ENSG_to_Name.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Jointure pour ajouter les noms de gènes
df_annotated <- left_join(df, Correspondance, by = c("Geneid" = "gene_id"))

# Réorganisation des colonnes (nom du gène juste après l’ID)
other_cols <- setdiff(colnames(df_annotated), c("Geneid", "gene_name"))
df_annotated <- df_annotated[, c("Geneid", "gene_name", other_cols)]

write.table(df_annotated, file = "all-counts_with_gene_names.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
############################################################################################################
#                               Étape 5 : Définition des gènes d’intérêt SLA                               #
############################################################################################################
genes_list <- c(
  "ALS2", "ANG", "ANO6", "APOE", "APP", "C21orf2", "C9orf72", "CAPRIN1", "CCNF", "CHCHD10",
  "CHMP2B", "CSF1R", "DAO", "DCTN1", "DPYSL3", "ELP3", "EPHA4", "ERBB4", "EWSR1", "FAS",
  "FIG4", "FUS", "GC", "GRN", "HNRNPA1", "HNRNPA2B1", "KIF5A", "LMNB1", "MAPT", "MATR3",
  "NEFH", "NEK1", "NR1H2", "NR1H3", "OPTN", "PFN1", "PRPH", "PSEN1", "PSEN2", "SETX",
  "SIGMAR1", "SOD1", "SPAST", "SPG11", "SPTLC1", "SQSTM1", "TAF15", "TARDBP", "TBK1", "TIA1",
  "TREM2", "TUBA4A", "UBQLN2", "UNC13A", "VAPB", "VCP"
)

############################################################################################################
#                         Étape 6 : Extraction des gènes d’intérêt dans la matrice annotée                 #
############################################################################################################
genes_filt <- df_annotated %>%
  filter(tolower(gene_name) %in% tolower(genes_list))

print(dim(genes_filt))
print("Controle du merge  ENSG -> Nom du gène")
print(head(genes_filt[, 1:2]))

# Sauvegarde finale
write.table(genes_filt, file = "counts_56genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

############################################################################################################
# Étape 7 : Définition des conditions expérimentales (si pas déjà définies)
############################################################################################################

# Extraction des noms d’échantillons (colonnes) dans counts_mat_filt
sample_names <- colnames(counts_mat_filt)

# Définition des conditions par nom d’échantillon, à adapter selon ton jeu de données
conditions <- ifelse(grepl("control", sample_names, ignore.case = TRUE), "control", "SLA")
conditions <- factor(conditions, levels = c("control", "SLA"))

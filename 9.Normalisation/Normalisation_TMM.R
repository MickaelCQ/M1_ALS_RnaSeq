############################################################################################################
# SCRIPT : RNA-seq SLA - Visualisation CPM (TMM) des 56 gènes d’intérêt
# Auteur : Mickael Coquerelle
# Date : 04/06/2025
# Description : Ce script lit une matrice de comptage RNA-seq ciblée,
#              calcule les CPM normalisés par la méthode TMM via edgeR,
#              puis génère des dotplots pour visualiser l’expression
#              normalisée des gènes d’intérêt par échantillon et condition.
############################################################################################################

############################################################################################################
#                           Chargement des librairies nécessaires
############################################################################################################
library(readr)     # Lecture de fichiers TSV
library(dplyr)     # Manipulation de données
library(tidyr)     # Transformation tidy (pivot)
library(tibble)    # Conversion des colonnes/rownames
library(edgeR)     # Calcul des CPM normalisés (TMM)
library(ggplot2)   # Visualisation

############################################################################################################
#                           Étape 1 : Lecture et mise en forme des données
############################################################################################################
# Chargement de la matrice de comptage des 56 gènes
counts <- read_tsv("~/Final_counts_56genes.tsv")

# Format long + séparation des métadonnées
counts_long <- counts %>%
  pivot_longer(cols = -c(Geneid, gene_name, Chr, Length),
               names_to = "full_id",
               values_to = "Counts") %>%
  separate(full_id, into = c("Run", "Sample", "Condition"), sep = "-", remove = TRUE)

# Agrégation des éventuels doublons (ex : isoformes fusionnées, ou comptages multiples)
counts_long_unique <- counts_long %>%
  group_by(gene_name, Sample) %>%
  summarise(Counts = sum(Counts, na.rm = TRUE)) %>%
  ungroup()

# Passage en format large (samples en colonnes)
counts_wide <- counts_long_unique %>%
  pivot_wider(names_from = Sample, values_from = Counts) %>%
  column_to_rownames(var = "gene_name")

# Conversion en matrice numérique
counts_mat <- as.matrix(counts_wide)
stopifnot(all(apply(counts_mat, 2, is.numeric)))

############################################################################################################
#                           Étape 2 : Normalisation TMM via edgeR
############################################################################################################
# Création de l'objet DGEList
dge <- DGEList(counts = counts_mat)

# Calcul des facteurs de normalisation TMM
dge <- calcNormFactors(dge, method = "TMM")

# Calcul des CPM normalisés
cpm_tmm <- cpm(dge, normalized.lib.sizes = TRUE)

############################################################################################################
#                       Étape 3 : Formatage tidy des données pour ggplot
############################################################################################################
# Passage en tidy data
cpm_tmm_df <- as.data.frame(cpm_tmm) %>%
  rownames_to_column("gene_name") %>%
  pivot_longer(cols = -gene_name, names_to = "Sample", values_to = "CPM_TMM")

# Ajout des informations Run et Condition
sample_info <- counts_long %>%
  distinct(Sample, Run, Condition)

cpm_tmm_long <- cpm_tmm_df %>%
  left_join(sample_info, by = "Sample")

############################################################################################################
#                 Étape 4 : Visualisation par groupes de gènes avec points colorés
############################################################################################################
# Liste des gènes
genes <- unique(cpm_tmm_long$gene_name)
gene_groups <- split(genes, ceiling(seq_along(genes) / 1)) # Groupes de 8 gènes

for (i in seq_along(gene_groups)) {
  gsub <- gene_groups[[i]]
  
  df_sub <- cpm_tmm_long %>%
    filter(gene_name %in% gsub)
  
  # Dotplot : un point par sample, coloré par Condition
  p <- ggplot(df_sub, aes(x = Run, y = CPM_TMM, color = Condition)) +
    geom_point(size = 2, alpha = 0.7, position = position_jitter(width = 0.25)) +
    facet_wrap(~ gene_name, scales = "free_y") +
    scale_color_brewer(palette = "Set2") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Expression CPM normalisée (TMM) - Groupe", i),
      x = "Runs",
      y = "CPM (normalisé TMM)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )
  
  print(p)
  readline(prompt = "Appuyez sur Entrée pour afficher le groupe suivant...")
}

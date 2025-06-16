############################################################################################################
# SCRIPT : RNA-seq SLA - Visualisation TPM des 56 gènes d’intérêt
# Auteur : Mickael Coquerelle
# Date : 04/06/2025
# Description : Ce script charge une matrice de comptage des 56 gènes SLA,
#              transforme les données en format long, calcule les TPM par échantillon,
#              puis génère des boxplots groupés par lots de 8 gènes pour visualiser
#              l’expression normalisée par Run et Condition.
############################################################################################################

############################################################################################################
#                           Chargement des librairies nécessaires
############################################################################################################
library(readr)     # Lecture de fichiers TSV
library(tidyverse) # Manipulation de données et visualisation (ggplot2, dplyr, tidyr)

############################################################################################################
#                           Étape 1 : Lecture des données et transformation en format long
############################################################################################################
# Chargement de la matrice de comptage des 56 gènes (fichier TSV avec colonnes de métadonnées et échantillons)
counts <- read_tsv("~/Final_counts_56genes.tsv")

# Passage en format "long" pour faciliter manipulation & visualisation
# - Toutes colonnes sauf Geneid, gene_name, Chr, Length sont pivotées en deux colonnes : full_id & Counts
# - Séparation de la colonne full_id en variables Run, Sample, Condition selon le séparateur "-"
counts_long <- counts %>%
  pivot_longer(
    cols = -c(Geneid, gene_name, Chr, Length),
    names_to = "full_id",
    values_to = "Counts"
  ) %>%
  separate(full_id, into = c("Run", "Sample", "Condition"), sep = "-", remove = TRUE)
############################################################################################################
#                                Étape 2 : Calcul des comptes normalisés TPM (Transcripts Per Million)
############################################################################################################
# Calcul du TPM pour chaque gène par échantillon (Sample)
# TPM = (Counts / Length) / somme (Counts/Length) * 1e6 pour normaliser la profondeur de séquençage et la longueur des gènes
counts_tpm <- counts_long %>%
  group_by(Sample) %>%
  mutate(
    rate = Counts / Length,
    tpm = rate / sum(rate) * 1e6
  ) %>%
  ungroup()

############################################################################################################
#                         Étape 3 : Préparation des groupes de gènes pour affichage par lots de 8
############################################################################################################
# Extraction de la liste unique des noms de gènes
genes <- unique(counts_tpm$gene_name)

# Découpage de la liste des gènes en groupes de 8 pour faciliter la lisibilité des plots
gene_groups <- split(genes, ceiling(seq_along(genes) / 4))

############################################################################################################
#                         Étape 4 : Boucle de visualisation des TPM par groupe de gènes
############################################################################################################
# Pour chaque groupe de gènes, création d’un boxplot montrant l’expression TPM selon Run et Condition (sla ou controle)
for(i in seq_along(gene_groups)) {
  gsub <- gene_groups[[i]]
  
  # Filtrage des données pour le groupe de gènes courant
  df_sub <- counts_tpm %>% filter(gene_name %in% gsub)
  
  # Création du graphique avec ggplot2
  p <- ggplot(df_sub, aes(x = Run, y = tpm, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +                   # Boxplot sans afficher les outliers
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +  # Ajout de points jitter pour visualiser la distribution
    facet_wrap(~ gene_name, scales = "free_y") +         # Facettes par gène avec échelles indépendantes
    scale_fill_brewer(palette = "Set2") +                # Palette de couleurs pour les conditions
    theme_minimal(base_size = 12) +                       # Thème minimal propre
    labs(
      title = paste("Expression TPM des gènes - Groupe", i),
      x = "Run",
      y = "TPM normalisé"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1), # Texte des axes X inclinés pour meilleure lisibilité
      strip.text = element_text(face = "bold"),          # Titre des facettes en gras
      legend.position = "top"                            # Légende en haut du graphique
    )
  print(p)
  
  # Pause dans la console interactive avant affichage du groupe suivant
  readline(prompt = "Appuyez sur Entrée pour afficher le groupe suivant...")
}

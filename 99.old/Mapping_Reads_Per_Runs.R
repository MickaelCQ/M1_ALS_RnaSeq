library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(scales)

#Clean:
rm(list = ls())

data <- read_csv("~/Stats_Log_merge.csv", show_col_types = FALSE)
head(data)

##########################################################
#  Mise au format long des datas pour simplifier le code #
########################################################## 

data_STAR <- data %>% pivot_longer(cols = c("STAR_Unique_reads", "STAR_Multi_reads", "STAR_No_map_reads"),
                                   names_to = "Alignement_STAR", values_to = "Reads_STAR")

data_CRAC <- data %>% pivot_longer(cols = c("CRAC_Unique_reads", "CRAC_Multi_reads", "CRAC_No_map_reads"),
                                   names_to = "Alignement_CRAC", values_to = "Reads_CRAC")

#################################
# Montage des labels de légende #
#################################
data_STAR$Alignement_STAR <- factor(data_STAR$Alignement_STAR,
                                    levels = c("STAR_Unique_reads", "STAR_Multi_reads", "STAR_No_map_reads"),
                                    labels = c("STAR Unique", "STAR Multiple", "STAR No map"))

data_CRAC$Alignement_CRAC <- factor(data_CRAC$Alignement_CRAC,
                                    levels = c("CRAC_Unique_reads","CRAC_Multi_reads","CRAC_No_map_reads"),
                                    labels = c("CRAC Unique","CRAC Multiple", "CRAC No map"))


###############################################################################
# Calcul des totaux pour définir les valeurs relatives par groupe : Run/Outil #
###############################################################################

# Sum_Reads : local, somme des lectures pour chaque combinaison Run et Alignement_STAR/CRAC
# Total     : somme des 'Sum_Reads' pour chaque Run.
# Sum_OUTIL : sous-dataframe du type -> [ Run |	Alignement_OUTIL |	Sum_Reads	 | Total	| Percent ]

# Pour chaque combo Run/Alignement on calcule la somme des reads (sous tableau)
Sum_STAR <- data_STAR %>% group_by(Run, Alignement_STAR) %>% summarise(Sum_Reads = sum(Reads_STAR), .groups = 'drop') %>% 
  # Somme total par Run cette fois pour calculer le pourcentage : 
  group_by(Run) %>% mutate(Total = sum(Sum_Reads), Percent = (Sum_Reads / Total) * 100) %>% ungroup()

Sum_CRAC <- data_CRAC %>% group_by(Run,Alignement_CRAC) %>% summarise(Sum_Reads = sum(Reads_CRAC), .groups = 'drop')%>%
  group_by(Run) %>% mutate(Total = sum(Sum_Reads), Percent =(Sum_Reads / Total) * 100) %>% ungroup()

#Sum_STAR
#Sum_CRAC

#############################################################
# Calcul de la profondeur globale pour chaque Run/Outil     #
#############################################################
Depth_STAR <- data %>% select(Run, STAR_Total_reads) %>% group_by(Run) %>% summarise(Depth = sum(STAR_Total_reads), .groups = 'drop')
Depth_CRAC <- data %>% select(Run, CRAC_Total_reads) %>% group_by(Run) %>% summarise(Depth = sum(CRAC_Total_reads), .groups = 'drop')

# Fusion des données de 'profondeur d'alignement' avec les données de somme
Join_STAR <- left_join(Sum_STAR, Depth_STAR, by = "Run")
Join_CRAC <- left_join(Sum_CRAC, Depth_CRAC, by ="Run")

# Ajout de la colonne Outil poour distinguer dans le fill alignement
Sum_STAR$Outil <- "STAR"
Sum_STAR <- Sum_STAR %>% rename(Alignement = Alignement_STAR)

Sum_CRAC$Outil <- "CRAC"
Sum_CRAC <- Sum_CRAC %>% rename(Alignement = Alignement_CRAC)

# Fusion
data_plot <- bind_rows(Sum_STAR, Sum_CRAC)

ggplot(data_plot, aes(x = Run, y = Sum_Reads, fill = Alignement, label = paste0(round(Percent, 1), "%"))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, aes(group = Outil)) +
  geom_text(aes(group = Outil), position = position_dodge(width = 0.8), vjust = 0.5, color = "white", size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("steelblue", "darkorange", "indianred","steelblue", "darkorange", "indianred")) +
  labs(
    title = "Évaluation comparative des alignements STAR et CRAC",
    subtitle = "Index : GRCh37 — STAR v2.7.8 & CRAC vX.X.X",
    x = "Run RNASeq (Année/Mois)",
    y = "Nombre de lectures alignées",
    fill = "Catégorie des lectures"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.y = element_text(angle = 45)
  ) +
  scale_y_continuous(labels = label_number(big.mark = ",", accuracy = 1))


# ggplot(Join_STAR, aes(x = Run, y = Sum_Reads, fill = Alignement, label = paste0(round(Percent, 1), "%"))) +
#   geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # passage en 'dodge' pour les barres côte-à-côte
#   geom_text(position = position_dodge(width = 0.7), color = "black", size = 3, vjust = -0.5) +
#   scale_fill_manual(values = c("steelblue", "darkorange", "indianred")) +
#   labs(title = "Distribution des lectures alignées par patient (GRCh37 avec STAR V2.7.8)",
#        x = "Run RNASeq (Année/Mois)", y = "Nombre de lectures", fill = "Catégorie",
#        subtitle = "Comparaison entre types d’alignement et variabilité inter-run") +
#   facet_wrap(~ Alignement, scales = "free_y") +  # Facette par type d’alignement
#   theme_minimal(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         plot.title = element_text(face = "bold", hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5),
#         strip.text = element_text(face = "bold")) +
#   scale_y_continuous(labels = label_number(big.mark = ",", accuracy = 1))
# Ajouter une courbe de tendance pour la 'profondeur d'alignement'


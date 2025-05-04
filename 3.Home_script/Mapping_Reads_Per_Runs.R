library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(scales)

rm(list = ls())

# Chargement des stats de Log BAM : 
data <- read_csv("~/Stats_Log_Final.csv", show_col_types = FALSE)

# Mise au format long pour se faciliter la vie dans les graphes ggplot2:
data_long <- data %>% pivot_longer(cols = c("Unique_reads", "Multi_reads", "No_map_reads"),
                      names_to = "Alignement", values_to = "Reads")

# Renommage en francais des sections d'alignement: 
data_long$Alignement <- factor(data_long$Alignement,
                               levels = c("Unique_reads", "Multi_reads", "No_map_reads"),
                               labels = c("Uniques", "Multiples", "Non-alignés"))

# Pour chaque combo Run/Alignement on calcule la somme des reads
data_sum <- data_long %>% group_by(Run, Alignement) %>%
  summarise(Sum_Reads = sum(Reads), .groups = 'drop') %>% group_by(Run) %>%
  mutate(Total = sum(Sum_Reads), Percent = (Sum_Reads / Total) * 100) %>% # Calcul du pourcentage global
  ungroup() # On dégroupe pour éviter les soucis dans les calculs suivants. 


# Calcul de la profondeur totale de lecture pour chaque Run
data_depth <- data %>%
  select(Run, Total_reads) %>%
  group_by(Run) %>%
  summarise(Depth = sum(Total_reads), .groups = 'drop')

# Fusion les données de 'profondeur d'alignement' avec les données de somme
data_combined <- left_join(data_sum, data_depth, by = "Run")

ggplot(data_combined, aes(x = Run, y = Sum_Reads, fill = Alignement, label = paste0(round(Percent, 1), "%"))) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(position = position_stack(vjust = 0.5), color = "white", size = 4, fontface = "bold") +
  scale_fill_manual(values = c("steelblue", "darkorange", "indianred")) +
  labs(
    title = "Distribution du mapping par run (GRCh37 avec STAR V2.7.8)",
    x = "Run RNASeq (Année/Mois)",
    y = "Nombre de lectures alignées",
    fill = "Catégorie"
  ) +  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.y = element_text(angle = 45)
  ) +
  scale_y_continuous(labels = label_number(big.mark = ",", accuracy = 1))

# ggplot(data_combined, aes(x = Run, y = Sum_Reads, fill = Alignement, label = paste0(round(Percent, 1), "%"))) +
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
  geom_smooth(data = data_combined, aes(x = Run, y = Depth), 
              method = "loess", color = "black", linetype ="solid" , size = 1, se = TRUE)


##############################################################
#           Chargement des librairies et nettoyage           #
##############################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(grid)
library(ggprism)
library(cowplot)
rm(list = ls())

##############################################################
#              Chargement et aperçu des données              #
##############################################################

data <- read_csv("~/Stats_Log_merge.csv", show_col_types = FALSE)
head(data)
##############################################################
#    Mise au format long des données STAR et CRAC (pivot)    #
##############################################################

data_STAR <- data %>% pivot_longer(
  cols = c("STAR_Unique_reads", "STAR_Multi_reads", "STAR_No_map_reads", "STAR_Chimeric_reads"),
  names_to = "Alignement_STAR", values_to = "Reads_STAR"
)

data_CRAC <- data %>% pivot_longer(
  cols = c("CRAC_Unique_reads", "CRAC_Multi_reads", "CRAC_No_map_reads", "CRAC_PE_Chimera_reads",),
  names_to = "Alignement_CRAC", values_to = "Reads_CRAC"
)

##############################################################
#         Definition des catégories avec des labels propres  #
##############################################################

data_STAR$Alignement_STAR <- factor(
  data_STAR$Alignement_STAR,
  levels = c("STAR_Unique_reads", "STAR_Multi_reads", "STAR_No_map_reads", "STAR_Chimeric_reads"),
  labels = c("STAR Unique", "STAR Multiple", "STAR No map", "STAR Chimeric")
)

data_CRAC$Alignement_CRAC <- factor(
  data_CRAC$Alignement_CRAC,
  levels = c("CRAC_Unique_reads", "CRAC_Multi_reads", "CRAC_No_map_reads", "CRAC_PE_Chimera_reads"),
  labels = c("CRAC Unique", "CRAC Multiple", "CRAC No map", "CRAC Chimeric")
)


##############################################################
#     Calcul du % de chaque catégorie par Run et Outil       #
##############################################################

Sum_STAR <- data_STAR %>%
  group_by(Run, Alignement_STAR) %>%
  summarise(Sum_Reads = sum(Reads_STAR), .groups = 'drop') %>%
  group_by(Run) %>%
  mutate(Total = sum(Sum_Reads), Percent = (Sum_Reads / Total) * 100) %>%
  ungroup()

Sum_CRAC <- data_CRAC %>%
  group_by(Run, Alignement_CRAC) %>%
  summarise(Sum_Reads = sum(Reads_CRAC)/2, .groups = 'drop') %>%
  group_by(Run) %>%
  mutate(Total = sum(Sum_Reads), Percent = (Sum_Reads / Total) * 100) %>%
  ungroup()

##############################################################
#         Ajout des profondeurs d’alignement globales        #
##############################################################

Depth_STAR <- data %>% select(Run, STAR_Total_reads) %>%
  group_by(Run) %>% summarise(Depth = sum(STAR_Total_reads), .groups = 'drop')

Depth_CRAC <- data %>% select(Run, CRAC_Total_reads) %>%
  group_by(Run) %>% summarise(Depth = sum(CRAC_Total_reads), .groups = 'drop')

Join_STAR <- left_join(Sum_STAR, Depth_STAR, by = "Run")
Join_CRAC <- left_join(Sum_CRAC, Depth_CRAC, by = "Run")

##############################################################
#        Fusion STAR/CRAC + Harmonisation des noms           #
##############################################################

Sum_STAR$Outil <- "STAR"
Sum_STAR <- Sum_STAR %>% rename(Alignement = Alignement_STAR)

Sum_CRAC$Outil <- "CRAC"
Sum_CRAC <- Sum_CRAC %>% rename(Alignement = Alignement_CRAC)

data_plot <- bind_rows(Sum_STAR, Sum_CRAC)

##############################################################
#       Reformatage pour couleur harmonisée & légende        #
##############################################################

data_plot <- data_plot %>%
  mutate(Category = case_when(
    grepl("Unique", Alignement)   ~ "Alignement unique",
    grepl("Multiple", Alignement) ~ "Alignement multiple",
    grepl("No map", Alignement)   ~ "No map",
    grepl("Chimeric", Alignement) ~ "Alignement chimeric"
  ))

data_plot$Category <- factor(data_plot$Category,
                             levels = c("Alignement unique", "Alignement multiple", "Alignement chimeric", "No map"))

couleurs_categories <- c(
  "Alignement unique"   = "steelblue",
  "Alignement multiple" = "orange",
  "Alignement chimeric" = "mediumorchid",
  "No map"              = "indianred"
)

##############################################################
#               Construction du graphique final              #
##############################################################

ggplot(data_plot, aes(x = Run, y = Sum_Reads, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, color = "gray30") +
  geom_text(
    aes(label = ifelse(Percent > 5, paste0(round(Percent, 1), ""), "")),
    position = position_stack(vjust = 0.5),
    size = 3, color = "black", fontface = "bold"
  ) +
  facet_wrap(~ Outil, nrow = 2) +
  scale_fill_manual(values = couleurs_categories, name = "Catégorie de lecture") +
  labs(
    title = "Répartition des lectures avec deix approches d'alignement STAR & CRAC",
    subtitle = "Index : GRCh37 — STAR v2.7.8 & CRAC v2.5.2",
    x = "Run RNASeq (Année/Mois)",
    y = "Répartitions des lectures exprimées en %",
    caption = ""
  ) +
  theme_bw() +
  
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.caption = element_text(hjust = 0.5, size = 12, face = "italic", color = "gray20"),
    axis.text.y = element_text(angle = 45, size = 10),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  scale_y_continuous(labels = label_number(big.mark = ",", accuracy = 1))

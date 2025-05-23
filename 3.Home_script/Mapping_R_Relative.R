###############################################################################
#                           CHARGEMENT DES LIBRAIRIES                         #
###############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
rm(list = ls())

###############################################################################
#                          IMPORT DES DONNÉES & COULEURS                      #
###############################################################################

data <- read_csv("~/Stats_Log_merge.csv", show_col_types = FALSE)

couleurs_categories <- c(
  "Unique"   = "#80b1d3",  # bleu clair
  "Multiple" = "#fdb462",  # orange clair
  "No_map"   = "#fb8072"   # rouge saumon
)

ordre_run <- unique(data$Run)

###############################################################################
#                         PRÉPARATION DES DONNÉES STAR                        #
###############################################################################

data_star <- data %>%
  group_by(Run) %>%
  summarise(
    Unique   = sum(STAR_Unique_reads),
    Multiple = sum(STAR_Multi_reads),
    No_map   = sum(STAR_No_map_reads + STAR_Chimeric_reads),
    .groups = "drop"
  ) %>%
  mutate(Outil = "STAR") %>%
  pivot_longer(cols = Unique:No_map, names_to = "Category", values_to = "Reads")

###############################################################################
#                         PRÉPARATION DES DONNÉES CRAC                        #
###############################################################################

data_crac <- data %>%
  group_by(Run) %>%
  summarise(
    Unique   = sum(CRAC_Unique_reads / 2),
    Multiple = sum(CRAC_Multi_reads / 2) + sum(CRAC_Dup_reads / 2),
    No_map   = sum(CRAC_No_map_reads / 2),
    .groups = "drop"
  ) %>%
  mutate(Outil = "CRAC") %>%
  pivot_longer(cols = Unique:No_map, names_to = "Category", values_to = "Reads")

###############################################################################
#                  FUSION, CALCUL DE POURCENTAGE & X_POS                      #
###############################################################################

data_plot <- bind_rows(data_star, data_crac) %>%
  group_by(Run, Outil) %>%
  mutate(Percent = 100 * Reads / sum(Reads)) %>%
  ungroup() %>%
  arrange(Run, Outil) %>%
  group_by(Run) %>%
  mutate(
    x_pos = as.numeric(factor(Run, levels = ordre_run)) * 2 +
      ifelse(Outil == "STAR", -0.42, 0.42)
  ) %>%
  ungroup()

# Réordonner les catégories pour l'affichage
data_plot <- data_plot %>%
  mutate(Category = factor(Category, levels = c("Unique", "No_map", "Multiple")))

###############################################################################
#                                GRAPHIQUE                                    #
###############################################################################

ggplot(data_plot, aes(x = x_pos, y = Percent, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", color = "gray30", width = 0.75) +
  scale_fill_manual(values = couleurs_categories) +
  scale_x_continuous(
    breaks = seq(2, by = 2, length.out = length(ordre_run)),
    labels = ordre_run
  ) +
  
  # Affichage des pourcentages (sauf < 2 %)
  geom_text(
    data = data_plot %>% filter(Category %in% c("Unique", "No_map")),
    aes(label = ifelse(Percent > 2, round(Percent, 1), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  
  # Affichage de l'outil sous les barres
  geom_text(
    data = data_plot %>% filter(Category == "Unique"),
    aes(y = -4, angle = 45, label = Outil),
    size = 3,
    color = "black"
  ) +
  
  labs(
    title    = "Étude comparative des alignements STAR et CRAC",
    subtitle = "Index : GRCh37 -- STAR v2.7.8 & CRAC v2.5.2",
    x        = "Run RNASeq (Année/Mois)",
    y        = "Pourcentage de lectures",
    caption  = "Répartitions des lectures exprimées en % (< 2 % non affichés)"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1.1),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(scale = 1),
    limits = c(-5, 100)
  )

###############################################################################
#                                EXPORT PDF                                   #
###############################################################################

ggsave("~/Figure_STAR_CRAC_R.pdf", width = 8, height = 5, dpi = 600)

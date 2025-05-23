library(dplyr)
library(tidyr)
library(ggplot2)
rm(list = ls())

##############################################################
#              Chargement et aperçu des données              #
##############################################################

data <- read_csv("~/Stats_Log_merge.csv", show_col_types = FALSE)

couleurs_categories <- c(
  "Unique"   = "#80b1d3",  # bleu clair
  "Multiple" = "#fdb462",  # orange clair
  "No_map"   = "#fb8072"   # rouge saumon
)

ordre_patient <- unique(data$Patient)

##############################################################
#        Préparation des données STAR par patient            #
##############################################################

data_star <- data %>%
  group_by(Patient) %>%
  summarise(
    Unique   = sum(STAR_Unique_reads),
    Multiple = sum(STAR_Multi_reads),
    No_map   = sum(STAR_No_map_reads + STAR_Chimeric_reads),
    .groups = "drop"
  ) %>%
  mutate(Outil = "STAR") %>%
  pivot_longer(cols = Unique:No_map, names_to = "Category", values_to = "Reads")

##############################################################
#        Préparation des données CRAC par patient            #
##############################################################

data_crac <- data %>%
  group_by(Patient) %>%
  summarise(
    Unique   = sum(CRAC_Unique_reads / 2),
    Multiple = sum(CRAC_Multi_reads / 2) + sum(CRAC_Dup_reads / 2),
    No_map   = sum(CRAC_No_map_reads / 2),
    .groups = "drop"
  ) %>%
  mutate(Outil = "CRAC") %>%
  pivot_longer(cols = Unique:No_map, names_to = "Category", values_to = "Reads")

##############################################################
#        Fusion et calcul des pourcentages par patient       #
##############################################################

data_plot <- bind_rows(data_star, data_crac) %>%
  group_by(Patient, Outil) %>%
  mutate(Percent = 100 * Reads / sum(Reads)) %>%
  ungroup()

data_plot <- data_plot %>%
  mutate(
    Patient = factor(Patient, levels = ordre_patient),
    Category = factor(Category, levels = c("Unique", "No_map", "Multiple"))
  )

##############################################################
#                   Génération du graphique                  #
##############################################################

ggplot(data_plot, aes(x = Patient, y = Percent, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", color = "gray30", width = 0.75) +
  scale_fill_manual(values = couleurs_categories) +
  scale_y_continuous(
    labels = scales::label_percent(scale = 1),
    limits = c(0, 101)
  ) +
  facet_grid(Outil ~ ., scales = "free_y") +
  geom_text(
    data = data_plot %>% filter(Category %in% c("Unique", "No_map")),
    aes(label = ifelse(Percent > 100, paste0(round(Percent, 1), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black",
    fontface = "bold"
  ) +
  labs(
    title = "Étude comparative des alignements STAR et CRAC par patient",
    subtitle = "Index : GRCh37 -- STAR v2.7.8 & CRAC v2.5.2",
    x = "Patient",
    y = "Pourcentage de lectures",
    fill = "Catégorie de lecture",
    caption = "Répartitions des lectures exprimées en %"
  ) +
  theme_minimal(base_size = 11, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9, face = "italic"),
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

##############################################################
#                     Sauvegarde du graphe                   #
##############################################################

ggsave("~/Figure_STAR_CRAC_P.pdf", width = 8, height = 5, dpi = 600)

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(scales)

# Lecture des données
data <- read_csv("~/Stats_Log_Final.csv", show_col_types = FALSE)

# Mise au format long
data_long <- data %>%
  pivot_longer(cols = c("Unique_reads", "Multi_reads", "No_map_reads"),
               names_to = "Alignement", values_to = "Reads")

# Renommage des niveaux d'alignement
data_long$Alignement <- factor(data_long$Alignement,
                               levels = c("Unique_reads", "Multi_reads", "No_map_reads"),
                               labels = c("Unique", "Multiple", "Non-aligné"))

# Calcul des pourcentages pour chaque patient
data_patient <- data_long %>%
  group_by(Patient, Alignement) %>%
  summarise(Sum_Reads = sum(Reads), .groups = "drop") %>%
  group_by(Patient) %>%
  mutate(Total = sum(Sum_Reads),
         Percent = (Sum_Reads / Total) * 100) %>%
  ungroup()

# Profondeur de lecture par patient
data_depth <- data %>% group_by(Patient) %>% summarise(Depth = sum(Total_reads), .groups = "drop")

# Ajout du champ Run pour chaque patient
patient_run <- data %>% select(Patient, Run) %>% distinct()

# Fusion des données
data_combined <- left_join(data_patient, data_depth, by = "Patient") %>%
  left_join(patient_run, by = "Patient")

# Ordre des patients dans le plot
data_combined$Patient <- factor(data_combined$Patient, 
                                levels = data_combined %>% arrange(Run, Patient) %>% pull(Patient) %>% unique())

# Obtenir les positions de début/fin par Run
run_ranges <- data_combined %>%
  distinct(Patient, Run) %>% mutate(Patient = factor(Patient, levels = levels(data_combined$Patient))) %>%
  arrange(Patient) %>% group_by(Run) %>%
  summarise(start = first(Patient), end = last(Patient)) %>%
  mutate(start_idx = as.numeric(start),
         end_idx = as.numeric(end),
         mid_idx = (start_idx + end_idx) / 2)

# Graphique principal
p <- ggplot(data_combined, aes(x = Patient, y = Sum_Reads, fill = Alignement)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_smooth(aes(x = as.numeric(Patient), y = Depth, group = 1), 
              method = "loess", se = TRUE, 
              color = "black", linetype = "solid", size = 0.5) +
  scale_fill_manual(values = c("steelblue", "darkorange", "indianred")) +
  scale_y_continuous(labels = label_number(big.mark = ",", accuracy = 1),
                     expand = expansion(mult = c(0.05, 0.01))) +
  labs(
    title = "Distribution du mapping/patient (GRCh37 avec STAR V2.7.8) ",
    subtitle = "avec courbe de tendance de la 'profondeur d'alignement (IC 95%) ",
    x = "Echantillons",
    y = "Nombre de reads ",
    fill = "Catégories"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.margin = margin(15, 15, 60, 20)
  )

# Ajout des annotations de groupe Run sous l'axe des patients
for (i in seq_len(nrow(run_ranges))) {
  p <- p +
    annotate("segment",
             x = run_ranges$start_idx[i], xend = run_ranges$end_idx[i],
             y = -max(data_combined$Sum_Reads) * 0.06,
             yend = -max(data_combined$Sum_Reads) * 0.06,
             size = 0.6, color = "black") +
    annotate("text",
             x = run_ranges$mid_idx[i],
             y = -max(data_combined$Sum_Reads) * 0.15,
             label = run_ranges$Run[i],
             size = 4, fontface = "italic")
}

# Affichage final
print(p)

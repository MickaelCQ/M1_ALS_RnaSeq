
# ===============================
# SCRIPT ANALYTIQUE : Heatmap corrélations - uniquement variables numériques
# ===============================

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

# 1. Chargement des données
fichier <- "Stats_Log_merge_with_deltas.csv"
df_raw <- read_csv(fichier, locale = locale(encoding = "ISO-8859-1"), guess_max = 5000, show_col_types = FALSE)

# 2. Sélection stricte des variables numériques parmi vars_exp et vars_resp
vars_exp <- c("Ville_Prescripteur", "Date_Prelevement", "Date_Recep", "Date_extraction", 
              "Concentration_ARN", "Purete_proteique", "Date_Lib", "Date_Lancement", 
              "Delta_Run_Prel", "Delta_Run_Recep", "Delta_Run_Ext", "Delta_Run_Lib", "STAR_Type")

vars_resp <- c("STAR_Total_reads", "STAR_Unique_reads", "STAR_Unique_pct", "STAR_Multi_reads", 
               "STAR_Multi_pct", "STAR_No_map_reads", "STAR_No_map_pct_sum",
               "CRAC_Total_reads", "CRAC_Unique_reads", "CRAC_Unique_pct", "CRAC_Multi_reads", 
               "CRAC_Multi_pct", "CRAC_No_map_reads", "CRAC_No_map_pct")

# Ne garder que les colonnes existantes dans df_raw
vars_exp <- intersect(vars_exp, colnames(df_raw))
vars_resp <- intersect(vars_resp, colnames(df_raw))

# Filtrer uniquement variables numériques explicatives (exclure dates et facteurs)
vars_exp_num <- vars_exp[sapply(df_raw[vars_exp], is.numeric)]

# Filtrer uniquement variables numériques réponses (exclure variables non numériques)
vars_resp_num <- vars_resp[sapply(df_raw[vars_resp], is.numeric)]

# 3. Fonction de calcul de corrélation
corr_metrics <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x_ok <- x[ok]
  y_ok <- y[ok]
  if(length(x_ok) < 10) return(c(cor = NA, R2 = NA, pval = NA))
  cor_val <- cor(x_ok, y_ok, method = "pearson")
  lmfit <- lm(y_ok ~ x_ok)
  sum_lm <- summary(lmfit)
  R2 <- sum_lm$r.squared
  pval <- coef(summary(lmfit))[2,4]
  c(cor = cor_val, R2 = R2, pval = pval)
}

# 4. Calcul des matrices de corrélation
cor_mat <- matrix(NA_real_, nrow = length(vars_resp_num), ncol = length(vars_exp_num),
                  dimnames = list(vars_resp_num, vars_exp_num))
R2_mat  <- cor_mat
pval_mat <- cor_mat

for (resp in vars_resp_num) {
  for (expv in vars_exp_num) {
    m <- corr_metrics(df_raw[[expv]], df_raw[[resp]])
    cor_mat[resp, expv] <- m["cor"]
    R2_mat[resp, expv] <- m["R2"]
    pval_mat[resp, expv] <- m["pval"]
  }
}

# 5. Construction dataframe tidy pour ggplot
df_corr <- as.data.frame(as.table(cor_mat))
colnames(df_corr) <- c("Response", "Explicative", "Correlation")
df_corr$R2 <- as.vector(R2_mat)
df_corr$pval <- as.vector(pval_mat)

df_corr <- df_corr %>%
  mutate(
    R2_label = ifelse(!is.na(R2), sprintf("R²=%.2f", R2), ""),
    signif = case_when(
      is.na(pval) ~ "NS",
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE        ~ "NS"
    ),
    label = ifelse(!is.na(Correlation),
                   paste0(sprintf("%.2f", Correlation), "\n", R2_label, "\n", signif),
                   "")
  )

# 6. Visualisation heatmap
heatmap_plot <- ggplot(df_corr, aes(x = Explicative, y = Response, fill = Correlation)) +
  geom_tile(color = "grey30") +
  geom_text(aes(label = label), size = 3, lineheight = 0.8, family = "mono") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0,
                       na.value = "grey90", limits = c(-1,1), name = "Corrélation (r)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "Heatmap des corrélations entre variables explicatives et réponses",
       subtitle = "Valeurs: corrélation r, R² et significativité (* p<0.05)",
       x = "Variables explicatives",
       y = "Variables réponses")

# 7. Sauvegarde & affichage
ggsave("Heatmap_Correlation_Explicatives_vs_Reponses_numeric_only.png", heatmap_plot, width = 14, height = 9, dpi = 300)
print(heatmap_plot)

# Export des résultats dans un CSV
write.csv(df_corr %>% select(Response, Explicative, Correlation, R2, pval, signif),
          file = "Correlation_Results_numeric_only.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

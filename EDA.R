### Plots für Einflussgrößen
source("preprocessing_final.R")
library(ggplot2)
library(ggtext)
library(knitr)
library(kableExtra)
library(reshape2)
library(heatmaply)
library(lsr)
library(effectsize)

# ----------------------------  Plots für Zielgrößen
# Verlaufsklassen mit Antibiotikagabe:
ggplot(data = lw_et_merged_with_all, aes(x = verlaufsklasse, fill = antibiotikagabe)) +
  ylab("Anzahl") + ggtitle("Absolute Häufigkeit der Zielgröße \"Verlaufsklasse\"") +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("palegreen4", "red3"),
                    name = "Antibiotikagabe") +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450)) + xlab("Verlaufsklasse") +
  theme(
    axis.text = element_text(
      size = 9,
      angle = 45,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )

# Antibiotika vs. kein Antibiotika
ggplot(data = lw_et_merged_with_all, aes(x = reorder(antibiotikagabe, antibiotikagabe, function(x)
  - length(x)))) +
  ylab("Anzahl") + ggtitle("Absolute Häufigkeit der Zielgröße \n \"Antibiotikagabe\"") +
  geom_bar(stat = "count", fill = c("red3", "palegreen4")) +
  ylim(c(0, 500)) +
  xlab("Antibiotikagabe") +
  theme(
    axis.text = element_text(
      size = 9,
      angle = 45,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title = element_text(size = 12, face = "bold")
  )


# Auffällig vs. Unauffällig:
tmp <- lw_et_merged_with_all
tmp$mit_auffaelligkeiten <-
  relevel(tmp$mit_auffaelligkeiten, ref = "ja")
tmp$mit_auffaelligkeiten <- as.factor(tmp$mit_auffaelligkeiten)

ggplot(data = tmp, aes(x = mit_auffaelligkeiten)) +
  ylab("Anzahl") + ggtitle("Absolute Häufigkeit der Zielgröße \n \"Auffälligkeiten auf dem Baum\"") +
  geom_bar(stat = "count", fill = c("plum3", "gold2")) +
  ylim(c(0, 600)) +
  xlab("Auffälligkeiten auf dem Baum") +
  theme(
    axis.text = element_text(
      size = 9,
      angle = 45,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title = element_text(size = 12, face = "bold")
  )

# ---------------------------- Korrelationen
# Korrelationen der metrischen Variablen mit metrischen Variablen nach Pearson:
metr_vars <-
  lw_et_merged %>% select(
    zellzahl,
    log_zellzahl,
    laktationsnummer,
    log_laktationsnummer,
    mean_milch,
    diff_milch,
    var_milch,
    last_milch,
    trend_zz
  )
cor_metr <- round(cor(metr_vars), 3)
melted_corr_mat <- melt(cor_metr)
# plotting the correlation heatmap
ggplot(melted_corr_mat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  coord_fixed() + geom_text(aes(Var2, Var1, label = value),
                            color = "black",
                            size = 4) +
  #label(x = "Einflussgrößen", y = "Einflussgrößen") +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Korrelations-\nkoeffizient nach \nBravais-Pearson"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 1,
    size = 12,
    hjust = 1
  )) +
  labs(y = "Einflussgrößen", x = "Einflussgrößen", title = "Heatmap für die Korrelationen zwischen den metrischen Einflussgrößen") +
  coord_fixed()

cor(as.numeric(lw_et_merged$zellzahl),
    as.numeric(lw_et_merged$diff_milch))
cor(as.numeric(lw_et_merged$zellzahl),
    as.numeric(lw_et_merged$mean_milch))
cor(as.numeric(lw_et_merged$zellzahl),
    as.numeric(lw_et_merged$var_milch))
cor(as.numeric(lw_et_merged$zellzahl),
    as.numeric(lw_et_merged$last_milch))
cor(as.numeric(lw_et_merged$zellzahl),
    as.numeric(lw_et_merged$trend_zz))

cor(as.numeric(lw_et_merged$diff_milch),
    as.numeric(lw_et_merged$mean_milch))
cor(as.numeric(lw_et_merged$diff_milch),
    as.numeric(lw_et_merged$var_milch))
cor(as.numeric(lw_et_merged$diff_milch),
    as.numeric(lw_et_merged$last_milch))
cor(as.numeric(lw_et_merged$diff_milch),
    as.numeric(lw_et_merged$trend_zz))

cor(as.numeric(lw_et_merged$mean_milch),
    as.numeric(lw_et_merged$trend_zz))
cor(as.numeric(lw_et_merged$mean_milch),
    as.numeric(lw_et_merged$var_milch))
cor(as.numeric(lw_et_merged$mean_milch),
    as.numeric(lw_et_merged$last_milch))

cor(as.numeric(lw_et_merged$var_milch),
    as.numeric(lw_et_merged$trend_zz))
cor(as.numeric(lw_et_merged$var_milch),
    as.numeric(lw_et_merged$last_milch))

cor(as.numeric(lw_et_merged$last_milch),
    as.numeric(lw_et_merged$trend_zz))

cor(log(as.numeric(lw_et_merged$zellzahl)), as.numeric(lw_et_merged$zellzahl))
cor(log(as.numeric(lw_et_merged$zellzahl)), as.numeric(lw_et_merged$diff_milch))
cor(log(as.numeric(lw_et_merged$zellzahl)), as.numeric(lw_et_merged$mean_milch))
cor(log(as.numeric(lw_et_merged$zellzahl)), as.numeric(lw_et_merged$var_milch))
cor(log(as.numeric(lw_et_merged$zellzahl)), as.numeric(lw_et_merged$last_milch))
cor(log(as.numeric(lw_et_merged$zellzahl)), as.numeric(lw_et_merged$trend_zz))

cor(log(as.numeric(lw_et_merged$laktationsnummer)), log(as.numeric(lw_et_merged$zellzahl)))
cor(log(as.numeric(lw_et_merged$laktationsnummer)), as.numeric(lw_et_merged$zellzahl))
cor(log(as.numeric(lw_et_merged$laktationsnummer)), as.numeric(lw_et_merged$diff_milch))
cor(log(as.numeric(lw_et_merged$laktationsnummer)), as.numeric(lw_et_merged$mean_milch))
cor(log(as.numeric(lw_et_merged$laktationsnummer)), as.numeric(lw_et_merged$var_milch))
cor(log(as.numeric(lw_et_merged$laktationsnummer)), as.numeric(lw_et_merged$last_milch))
cor(log(as.numeric(lw_et_merged$laktationsnummer)), as.numeric(lw_et_merged$trend_zz))

cor(as.numeric(lw_et_merged$laktationsnummer), log(as.numeric(lw_et_merged$zellzahl)))
cor(as.numeric(lw_et_merged$laktationsnummer),
    as.numeric(lw_et_merged$zellzahl))
cor(as.numeric(lw_et_merged$laktationsnummer),
    as.numeric(lw_et_merged$diff_milch))
cor(as.numeric(lw_et_merged$laktationsnummer),
    as.numeric(lw_et_merged$mean_milch))
cor(as.numeric(lw_et_merged$laktationsnummer),
    as.numeric(lw_et_merged$var_milch))
cor(as.numeric(lw_et_merged$laktationsnummer),
    as.numeric(lw_et_merged$last_milch))
cor(as.numeric(lw_et_merged$laktationsnummer),
    as.numeric(lw_et_merged$trend_zz))
cor(as.numeric(lw_et_merged$laktationsnummer), log(as.numeric(lw_et_merged$laktationsnummer)))

# Korrelationen der metrischen Variablen mit ordinalen Variablen (kategorie) nach Spearman:
cor_spearman <-
  cor(metr_vars,
      as.numeric(lw_et_merged$kategorie),
      method = "spearman",
      use = "complete.obs")

cor(
  as.numeric(lw_et_merged$zellzahl),
  as.numeric(lw_et_merged$kategorie),
  method = "spearman",
  use = "complete.obs"
)

cor(
  as.numeric(lw_et_merged$diff_milch),
  as.numeric(lw_et_merged$kategorie),
  method = "spearman",
  use = "complete.obs"
)

cor(
  as.numeric(lw_et_merged$mean_milch),
  as.numeric(lw_et_merged$kategorie),
  method = "spearman",
  use = "complete.obs"
)

cor(
  as.numeric(lw_et_merged$var_milch),
  as.numeric(lw_et_merged$kategorie),
  method = "spearman",
  use = "complete.obs"
)

cor(
  as.numeric(lw_et_merged$last_milch),
  as.numeric(lw_et_merged$kategorie),
  method = "spearman",
  use = "complete.obs"
)

cor(
  as.numeric(lw_et_merged$trend_zz),
  as.numeric(lw_et_merged$kategorie),
  method = "spearman",
  use = "complete.obs"
)

cor(
  log(as.numeric(lw_et_merged$zellzahl)),
  as.numeric(lw_et_merged$kategorie),
  method = "spearman",
  use = "complete.obs"
)

cor(
  log(as.numeric(lw_et_merged$laktationsnummer)),
  as.numeric(lw_et_merged$kategorie),
  method = "spearman",
  use = "complete.obs"
)

cor(
  as.numeric(lw_et_merged$laktationsnummer),
  as.numeric(lw_et_merged$kategorie),
  method = "spearman",
  use = "complete.obs"
)


# Korrelationen der nominalen Variablen (jahreszeit) mit metrischen Variablen nach Eta-Quadrat:
etaSquared(x = lm((
  lw_et_merged$zellzahl ~ lw_et_merged$jahreszeit
)))
etaSquared(x = lm((
  log(lw_et_merged$zellzahl) ~ lw_et_merged$jahreszeit
)))
etaSquared(x = lm((
  lw_et_merged$diff_milch ~ lw_et_merged$jahreszeit
)))
etaSquared(x = lm((
  lw_et_merged$mean_milch ~ lw_et_merged$jahreszeit
)))
etaSquared(x = lm((
  lw_et_merged$var_milch ~ lw_et_merged$jahreszeit
)))
etaSquared(x = lm((
  lw_et_merged$last_milch ~ lw_et_merged$jahreszeit
)))
etaSquared(x = lm((
  lw_et_merged$trend_zz ~ lw_et_merged$jahreszeit
)))
etaSquared(x = lm((
  log(lw_et_merged$laktationsnummer) ~ lw_et_merged$jahreszeit
)))
etaSquared(x = lm((
  lw_et_merged$laktationsnummer ~ lw_et_merged$jahreszeit
)))

eta_squared(lm(lw_et_merged$zellzahl ~ lw_et_merged$jahreszeit), partial = FALSE)
eta_squared(lm(log(lw_et_merged$zellzahl) ~ lw_et_merged$jahreszeit), partial = FALSE)
eta_squared(lm(lw_et_merged$diff_milch ~ lw_et_merged$jahreszeit),
            partial = FALSE)
eta_squared(lm(lw_et_merged$mean_milch ~ lw_et_merged$jahreszeit),
            partial = FALSE)
eta_squared(lm(lw_et_merged$var_milch ~ lw_et_merged$jahreszeit),
            partial = FALSE)
eta_squared(lm(lw_et_merged$last_milch ~ lw_et_merged$jahreszeit),
            partial = FALSE)
eta_squared(lm(lw_et_merged$trend_zz ~ lw_et_merged$jahreszeit), partial = FALSE)
eta_squared(lm(log(lw_et_merged$laktationsnummer) ~ lw_et_merged$jahreszeit), partial = FALSE)
eta_squared(lm(lw_et_merged$laktationsnummer ~ lw_et_merged$jahreszeit),
            partial = FALSE)


# Korrelationen der dichotomen Variablen (Schalmtest) mit metrischen Variablen nach Punktbiseraler Korrelation:
cor.test(as.numeric(lw_et_merged$schalmtest),
         as.numeric(lw_et_merged$zellzahl))
cor.test(as.numeric(lw_et_merged$schalmtest), log(as.numeric(lw_et_merged$zellzahl)))
cor.test(as.numeric(lw_et_merged$schalmtest),
         as.numeric(lw_et_merged$diff_milch))
cor.test(as.numeric(lw_et_merged$schalmtest),
         as.numeric(lw_et_merged$mean_milch))
cor.test(as.numeric(lw_et_merged$schalmtest),
         as.numeric(lw_et_merged$var_milch))
cor.test(as.numeric(lw_et_merged$schalmtest),
         as.numeric(lw_et_merged$last_milch))
cor.test(as.numeric(lw_et_merged$schalmtest),
         as.numeric(lw_et_merged$trend_zz))
cor.test(as.numeric(lw_et_merged$schalmtest), log(as.numeric(lw_et_merged$laktationsnummer))) # 0.03126443
cor.test(as.numeric(lw_et_merged$schalmtest),
         as.numeric(lw_et_merged$laktationsnummer))

source("preprocessing_final.R")
library(nnet)
library(brms)
library(lme4)
library(pROC)
library(ROCR)
library(caret)
library(tidyverse)
library(mclogit)
library(VGAM)
library(EffectStars2)
library(mboost)
library(glmnet)
library(splines)
library(boot)
library(gamclass)
library(gridExtra)
set.seed(123123)
### Hilfsfunktionen:
# Funktion, die Spezifikationen mit hochkorrelierten Variablen aussortiert
# Input:
# df = Dataframe mit Spezifikation (1 Spezfiktaion == 1 Spalte)
remove_combn <- function(df) {
  df <-
    df[, which(apply(df, 2, function(x)
      ! (("diff_milch" %in% x) & ("var_milch" %in% x)
      )))]
  df <-
    df[, which(apply(df, 2, function(x)
      ! (("zellzahl" %in% x) & ("log_zellzahl" %in% x)
      )))]
  df <-
    df[, which(apply(df, 2, function(x)
      ! (("laktationsnummer" %in% x) &
           ("log_laktationsnummer" %in% x)
      )))]
  df
}

# Funktion für die Kreuzvalidierung
# Input:
# k = Anzahl der Folds
# formula = Modellgleichung
# data = Daten
# target = Zielgröße
cross_validation_bin <- function(k, formula, data, target) {
  set.seed(123123)
  cv_auc <- rep(0, k)
  # create fold-indices
  folds <- cut(seq(1, nrow(data)), breaks = k, labels = FALSE)
  data$fold <- folds
  for (i in 1:k) {
    train <- filter(data, fold != i)
    test <- filter(data, fold == i)
    # Modell auf Trainingsdatensatz schaetzen
    fitted.mod <- glm(formula, data = train, family = "binomial")
    # Wahrscheinichkeit für Testdaten vorhersagen
    pred_probs <-
      predict(fitted.mod, type = "response", newdata = test)
    cv_auc[i] <-
      roc(test$antibiotikagabe, as.vector(pred_probs))$auc[1]

  }
  c(Reduce(paste, deparse(formula)), mean(cv_auc, na.rm = T))
}

# ---------------------------- Alle Formeln mit allen möglichen Variablen-Kombinationen erstellen
vars <-
  c(
    "mean_milch",
    "diff_milch",
    "last_milch",
    "trend_zz",
    "zellzahl",
    "log_zellzahl",
    "log_laktationsnummer",
    "jahreszeit",
    "var_milch",
    "kategorie",
    "schalmtest",
    "laktationsnummer"
  )

c9 <- combn(vars, m = 9)
c8 <- combn(vars, m = 8)
c7 <- combn(vars, m = 7)
c6 <- combn(vars, m = 6)
c5 <- combn(vars, m = 5)
c4 <- combn(vars, m = 4)
c3 <- combn(vars, m = 3)
c2 <- combn(vars, m = 2)
c1 <- combn(vars, m = 1)
formulas9 <-
  apply(remove_combn(c9), 2, function(x)
    formula(paste(
      "antibiotikagabe", paste(x, collapse = "+"), sep = "~"
    )))
formulas8 <-
  apply(remove_combn(c8), 2, function(x)
    formula(paste(
      "antibiotikagabe", paste(x, collapse = "+"), sep = "~"
    )))
formulas7 <-
  apply(remove_combn(c7), 2, function(x)
    formula(paste(
      "antibiotikagabe", paste(x, collapse = "+"), sep = "~"
    )))
formulas6 <-
  apply(remove_combn(c6), 2, function(x)
    formula(paste(
      "antibiotikagabe", paste(x, collapse = "+"), sep = "~"
    )))
formulas5 <-
  apply(remove_combn(c5), 2, function(x)
    formula(paste(
      "antibiotikagabe", paste(x, collapse = "+"), sep = "~"
    )))
formulas4 <-
  apply(remove_combn(c4), 2, function(x)
    formula(paste(
      "antibiotikagabe", paste(x, collapse = "+"), sep = "~"
    )))
formulas3 <-
  apply(remove_combn(c3), 2, function(x)
    formula(paste(
      "antibiotikagabe", paste(x, collapse = "+"), sep = "~"
    )))
formulas2 <-
  apply(remove_combn(c2), 2, function(x)
    formula(paste(
      "antibiotikagabe", paste(x, collapse = "+"), sep = "~"
    )))
formulas1 <-
  apply(c1, 2, function(x)
    formula(paste(
      "antibiotikagabe", paste(x, collapse = "+"), sep = "~"
    )))

# ---------------------------- Dataframes für CV-Ergebnisse für alle Kombinationen erstellen
l9 <- length(formulas9)
df_performance9_bin <-
  data.frame(
    "nmb of features" = rep(9, l9),
    "formula" = rep(0, l9),
    "cv_auc" = rep(0, l9)
  )
l8 <- length(formulas8)
df_performance8_bin <-
  data.frame(
    "nmb of features" = rep(8, l8),
    "formula" = rep(0, l8),
    "cv_auc" = rep(0, l8)
  )
l7 <- length(formulas7)
df_performance7_bin <-
  data.frame(
    "nmb of features" = rep(7, l7),
    "formula" = rep(0, l7),
    "cv_auc" = rep(0, l7)
  )
l6 <- length(formulas6)
df_performance6_bin <-
  data.frame(
    "nmb of features" = rep(6, l6),
    "formula" = rep(0, l6),
    "cv_auc" = rep(0, l6)
  )
l5 <- length(formulas5)
df_performance5_bin <-
  data.frame(
    "nmb of features" = rep(5, l5),
    "formula" = rep(0, l5),
    "cv_auc" = rep(0, l5)
  )
l4 <- length(formulas4)
df_performance4_bin <-
  data.frame(
    "nmb of features" = rep(4, l4),
    "formula" = rep(0, l4),
    "cv_auc" = rep(0, l4)
  )
l3 <- length(formulas3)
df_performance3_bin <-
  data.frame(
    "nmb of features" = rep(3, l3),
    "formula" = rep(0, l3),
    "cv_auc" = rep(0, l3)
  )
l2 <- length(formulas2)
df_performance2_bin <-
  data.frame(
    "nmb of features" = rep(2, l2),
    "formula" = rep(0, l2),
    "cv_auc" = rep(0, l2)
  )
l1 <- length(formulas1)
df_performance1_bin <-
  data.frame(
    "nmb of features" = rep(1, l1),
    "formula" = rep(0, l1),
    "cv_auc" = rep(0, l1)
  )


# ---------------------------- CV für alle Kombinationen durchführen
for (i in seq(formulas9)) {
  cv <-
    cross_validation_bin(10, formulas9[[i]], lw_et_merged, "antibiotikagabe")
  df_performance9_bin[i, c(2:3)] <- cv
}

for (i in seq(formulas8)) {
  cv <-
    cross_validation_bin(10, formulas8[[i]], lw_et_merged, "antibiotikagabe")
  df_performance8_bin[i, c(2:3)] <- cv
}

for (i in seq(formulas7)) {
  cv <-
    cross_validation_bin(10, formulas7[[i]], lw_et_merged, "antibiotikagabe")
  df_performance7_bin[i, c(2:3)] <- cv
}

for (i in seq(formulas6)) {
  cv <-
    cross_validation_bin(10, formulas6[[i]], lw_et_merged, "antibiotikagabe")
  df_performance6_bin[i, c(2:3)] <- cv
}

for (i in seq(formulas5)) {
  cv <-
    cross_validation_bin(10, formulas5[[i]], lw_et_merged, "antibiotikagabe")
  df_performance5_bin[i, c(2:3)] <- cv
}

for (i in seq(formulas4)) {
  cv <-
    cross_validation_bin(10, formulas4[[i]], lw_et_merged, "antibiotikagabe")
  df_performance4_bin[i, c(2:3)] <- cv
}

for (i in seq(formulas3)) {
  cv <-
    cross_validation_bin(10, formulas3[[i]], lw_et_merged, "antibiotikagabe")
  df_performance3_bin[i, c(2:3)] <- cv
}

for (i in seq(formulas2)) {
  cv <-
    cross_validation_bin(10, formulas2[[i]], lw_et_merged, "antibiotikagabe")
  df_performance2_bin[i, c(2:3)] <- cv
}

for (i in seq(formulas1)) {
  cv <-
    cross_validation_bin(10, formulas1[[i]], lw_et_merged, "antibiotikagabe")
  df_performance1_bin[i, c(2:3)] <- cv
}

# ---------------------------- Ergebnisse der CV
df_overview_bin <-
  rbind(
    df_performance9_bin,
    df_performance8_bin,
    df_performance7_bin,
    df_performance6_bin,
    df_performance5_bin,
    df_performance4_bin,
    df_performance3_bin,
    df_performance2_bin,
    df_performance1_bin
  )
df_overview_ordered_bin_auc <-
  df_overview_bin[order(df_overview_bin$cv_auc, decreasing = TRUE), ]

# ---------------------------- Stepwise AIC zur explorative Identifikation von Interaktionen
formula_bin_max_auc <-
  formula(df_overview_ordered_bin_auc[1, "formula"])
string_bin_upper_interaction <-
  paste0("(", paste(all.vars(formula_bin_max_auc[-2]), collapse = " + "), ")^2")
formula_bin_upper_interaction <-
  formula(paste("antibiotikagabe", string_bin_upper_interaction , sep = " ~ "))


mod_interaktion_bin_upper <- glm(formula_bin_upper_interaction,
                                 data = lw_et_merged,
                                 family = "binomial")
mod_interaktion_bin_lower <- glm(formula_bin_max_auc,
                                 lw_et_merged, family = "binomial")
interaktion_bin_sel <-
  step(
    mod_interaktion_bin_upper,
    direction = "both",
    trace = 0,
    scope = list(upper = mod_interaktion_bin_upper, lower = mod_interaktion_bin_lower)
  )

df_performance_bin_sel <-
  data.frame(
    "nmb of features" = c(6),
    "formula" = c(8),
    "cv_auc" = c(8)
  )
cv <-
  cross_validation_bin(10,
                       formula(interaktion_bin_sel$call),
                       lw_et_merged,
                       "antibiotikagabe")
df_performance_bin_sel[1, c(2:3)] <- cv

# ---------------------------- Finale Ergebnisse
df_overview_bin <-
  rbind(
    df_performance9_bin,
    df_performance8_bin,
    df_performance7_bin,
    df_performance6_bin,
    df_performance5_bin,
    df_performance4_bin,
    df_performance_bin_sel
  )
df_overview_ordered_bin_auc <-
  df_overview_bin[order(df_overview_bin$cv_auc, decreasing = TRUE), ]

pdf("Overview_bin.pdf", height = 30, width = 40)
grid.table(df_overview_ordered_bin_auc[c(1:100),])
dev.off()

save(df_overview_ordered_bin_auc, file = "overview_bin_final.RData")

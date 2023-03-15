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
library(effects)
library(gamclass)
library(gridExtra)
library(tibble)
library(sjPlot)
set.seed(123123)

#Hilfsfunktionen
lw_et_merged_af <- lw_et_merged_with_all[lw_et_merged_with_all$verlaufsklasse == 6, ]
cross_validation_bin_final <- function(k, formula, data, target){
  set.seed(123123)
  cv_f1 <- rep(0, k)
  cv_accuracy <- rep(0, k)
  cv_accuracy_balanced <- rep(0, k)
  cv_sens <- rep(0, k)
  cv_spec <- rep(0, k)
  cv_auc <- rep(0,k)
  cv_auc_global <- rep(0,k)
  warnings_auc <- rep(0,k)
  data <- data[sample(nrow(data)),]

  # create fold-indices
  folds <- cut(seq(1, nrow(data)), breaks = k, labels = FALSE)
  data$fold <- folds
  for (i in 1:k){
    # Train und Testdatensatz definieren
    train <- filter(data, fold != i)
    test <- filter(data, fold == i)
    # Modell auf Trainingsdatensatz schaetzen
    fitted.mod <- glm(formula, data = train, family = "binomial")
    # Wahrscheinichkeit für Testdaten vorhersagen
    pred_probs <- predict(fitted.mod, type = "response", newdata=test)
    pred_class <- ifelse(pred_probs > 0.5, "ja", "nein")
    # Performancemetriken berechnen
    cm <- confusionMatrix(as.factor(pred_class),
                          as.factor(test$antibiotikagabe))
    cv_accuracy[i] <- cm$overall[[1]]
    cv_accuracy_balanced[i] <- cm$byClass[[11]]
    cv_f1[i] <- cm$byClass[[7]]
    cv_sens[i] <- cm$byClass[[1]]
    cv_spec[i] <- cm$byClass[[2]]
    cv_auc[i] <- roc(test$antibiotikagabe, as.vector(pred_probs))$auc[1]

    # Globale Schätzer
    test_global_antibiotikagabe <- c(test$antibiotikagabe, lw_et_merged_af$antibiotikagabe)
    pred_probs_global <- c(pred_probs, rep(1, nrow(lw_et_merged_af)))
    pred_class_global <- c(pred_class, rep("ja", nrow(lw_et_merged_af)))
    cv_auc_global[i] <- roc(test_global_antibiotikagabe, pred_probs_global)$auc[1]
    # cm_global <- confusionMatrix(as.factor(pred_class_global),
    #                              as.factor(test_global_antibiotikagabe))
    # print(cm_global)
    # print(cm)
    # print(cm_global$byClass[[11]])
  }
  c(Reduce(paste, deparse(formula)), mean(cv_auc, na.rm = T), mean(cv_auc_global, na.rm = T),
    mean(cv_f1, na.rm = T), mean(cv_accuracy_balanced, na.rm = T), mean(cv_accuracy, NA.rm = T),
    mean(cv_sens, na.rm = T), mean(cv_spec, na.rm = T))
}

# ---------------------------- Modell und Formel laden
#source("modellwahl_multinomial")
load("overview_bin_final.RData")
final_model_bin_formula <- formula(df_overview_ordered_bin_auc[1, "formula"])
#formula_plot <- formula(antibiotikagabe ~ diff_milch + log_zellzahl + last_milch + log_laktationsnummer +
#schalmtest + log_zellzahl*last_milch)
# ---------------------------- Bestes Modell schätzen
final_model_bin <-
  glm(
    final_model_bin_formula,
    data = lw_et_merged, family = "binomial"(link = "logit"))

# ---------------------------- Performance schätzen
performance_final_model_bin <- data.frame("formula" = c(8), "cv_auc" = c(8) , "cv_auc_global" = c(8),
                                            "cv_f1"= c(8), "cv_accuarcy_balanced"= c(8),
                                            "cv_accuarcy"= c(8), "cv_sensitivitaet"= c(8), "cv_spezifitaet"= c(8))
performance_final_model_bin[1,] <- cross_validation_bin_final(10, tmp, lw_et_merged, "antibiotikagabe")
# ---------------------------- Forest- und Effektplot erstellen
forestplot <- plot_model(final_model_bin, value.offset = .3, axis.title = "Exponierte Koeffizientenschätzer",
           sort.est = T, axis.lim = c(0.01, 20), show.values = TRUE,
           title = "Forest-Plot der Koeffizientenschätzer der Einflussgrößen")#, show.intercept = TRUE)
forestplot
effectplot <- plot_model(final_model_bin, type = "int",
                         mdrt.values = "meansd", terms=c("log_zellzahl", "last_milch"),
                         title = "Effektplot der Interaktion zwischen log. Zellzahl und Abmelkungsmilchmenge",
                         axis.title = "Vorhergesagte Wahrscheinlichkeit für Antibiotikagabe \"ja\"")
effectplot

#knitr::kable(summary(final_model_bin)$coefficients)
#plot_model(final_model_bin, type = "int", mdrt.values = "meansd", terms=c("log_zellzahl", "last_milch [all]"))
#summary(final_model_bin)$coefficients

#library(stargazer)
#stargazer::stargazer(final_model_bin, title = "Modelloutput binäre logsitische Regression")
#stargazer::stargazer(final_model_multi, title = "Modelloutput multinomiale logistische Regression (ohne Betriebseffekt)")
#stargazer::stargazer(final_model_multi_rnd_eff$coefficients)
#final_model_multi_rnd_eff

library(sjstats)
library(VGAM)
library(EffectStars)
library(gt)
library(webshot2)
library(tidyverse)
library(nnet)
library(caret)
library(pROC)
library(mclogit)
source("preprocessing_final.R")

# ---------------------------- Hilfsfunktionen für finale Performanceschätzung
# Anzahl an Tieren aus Verlaufsklasse 6 damit der Anteil im globalen Testset stimmt
lw_et_merged_af <- lw_et_merged_with_all[lw_et_merged_with_all$verlaufsklasse == 6, ]
lw_et_merged_af <- lw_et_merged_af[c(1:37), ]
l_af <- 37
pred_af <- data.frame("1" = rep(0,l_af), "2" = rep(0,l_af), "3" = rep(0,l_af) ,
                      "4"= rep(0,l_af), "5"= rep(0,l_af), "6"= rep(1,l_af), prob_antibiotikagabe = rep(1,l_af), antibiotikagabe = rep(1,l_af))
colnames(pred_af) <- c("1","2", "3", "4","5","6", "prob_antibiotikagabe", "antibiotikagabe")

### Performance für multinomiales Modell ohne zufälligem Betriebseffekt
cross_validation_multinom_final <- function(k, formula, data, target){
  set.seed(123123)
  cv_f1 <- rep(0, k)
  cv_accuarcy <- rep(0, k)
  cv_accuarcy_balanced <- rep(0, k)
  cv_macro_f1 <- rep(0, k)
  cv_kappa <- rep(0, k)
  cv_auc <- rep(0,k)
  cv_auc_global <- rep(0,k)
  cv_auc_binary <- rep(0,k)
  cv_auc_binary_global <- rep(0,k)
  warnings_auc <- rep(0,k)
  data<-data[sample(nrow(data)),]
  cv_accuracy_bin <- rep(0,k)
  cv_accuracy_balanced_bin <- rep(0,k)
  cv_f1_bin <- rep(0,k)
  cv_sens_bin <- rep(0,k)
  cv_spec_bin <- rep(0,k)

  # create fold-indices
  folds <- cut(seq(1,nrow(data)),breaks=k,labels=FALSE)
  data$fold <- folds
  for (i in 1:k){
    # Train und Testdatensatz definieren
    train <- filter(data, fold != i)
    test <- filter(data, fold == i)
    # Modell auf Trainingsdatensatz schaetzen
    fitted.mod <- multinom(formula, data = train)
    # Wahrscheinichkeit für Testdaten vorhersagen
    pred_probs <- predict(fitted.mod, type = "prob", newdata=test)
    pred_class <- predict(fitted.mod, type = "class", newdata=test)
    # Performancemetriken berechnen
    cm <- confusionMatrix(pred_class, test[,target])
    cv_accuarcy[i] <- cm$overall[[1]]
    cv_accuarcy_balanced[i] <- sum(cm$byClass[,6], na.rm = T)/5
    cv_macro_f1[i] <- 2 * (sum(cm$byClass[, 5], na.rm = T)/5 * sum(cm$byClass[, 6], na.rm = T)/5)/(5/sum(cm$byClass[, 5], na.rm = T) + 5/sum(cm$byClass[,6], NA.rm = T))
    cv_kappa[i] <- cm$overall[[2]]
    cv_auc[i] <- multiclass.roc(test[,target], pred_probs)$auc[1]

    ### Wsl. fuer binaere Zielgroeße antibiotikagabe hinzufuegen
    pred_probs <- as.data.frame(pred_probs)
    pred_probs$prob_antibiotikagabe <- rowSums(pred_probs[3:5])
    pred_probs$antibiotikagabe <- ifelse(pred_class %in% c(1,2), "nein", "ja")
    cv_auc_binary[i] <- auc(roc(test$antibiotikagabe, pred_probs$prob_antibiotikagabe))
    cm_binaer <- confusionMatrix(as.factor(pred_probs$antibiotikagabe), test$antibiotikagabe)
    cv_accuracy_bin[i] <- cm_binaer$overall[[1]]
    cv_accuracy_balanced_bin[i] <- cm_binaer$byClass[[11]]
    cv_f1_bin[i] <- cm_binaer$byClass[[7]]
    cv_sens_bin[i] <- cm_binaer$byClass[[1]]
    cv_spec_bin[i] <- cm_binaer$byClass[[2]]

    # Globale Schätzer
    test_global <- c(test[,target], lw_et_merged_af$verlaufsklasse)
    test_global_antibiotikagabe <- c(test$antibiotikagabe, lw_et_merged_af$antibiotikagabe)
    pred_probs_global <- add_column(pred_probs, "6" = rep(0, nrow(pred_probs)), .after = "5")
    pred_probs_global <- rbind(pred_probs_global, pred_af)
    cv_auc_global[i] <- multiclass.roc(test_global, pred_probs_global[,  c("1","2", "3", "4","5", "6")])$auc[1]
    cv_auc_binary_global[i] <- auc(roc(test_global_antibiotikagabe, pred_probs_global$prob_antibiotikagabe))
  }
  c(Reduce(paste, deparse(formula)), mean(cv_auc, na.rm = T), mean(cv_auc_global, na.rm = T),
    mean(cv_auc_binary, na.rm = T), mean(cv_auc_binary_global, na.rm = T), mean(cv_macro_f1, na.rm = T),
    mean(cv_accuarcy_balanced, na.rm = T), mean(cv_accuarcy, NA.rm = T), mean(cv_kappa, NA.rm = T),
    mean(cv_accuracy_bin, na.rm = T), mean(cv_accuracy_balanced_bin, na.rm = T),
    mean(cv_f1_bin, na.rm = T), mean(cv_sens_bin, na.rm = T), mean(cv_spec_bin, na.rm = T))
}

### Performance für multinomiales Modell mit zufälligem Betriebseffekt
cross_validation_multinom_rnd_effect_final <- function(k, formula, data, target){
  set.seed(123123)
  cv_f1 <- rep(0, k)
  cv_accuarcy <- rep(0, k)
  cv_accuarcy_balanced <- rep(0, k)
  cv_macro_f1 <- rep(0, k)
  cv_kappa <- rep(0, k)
  cv_auc <- rep(0,k)
  cv_auc <- rep(0,k)
  cv_auc_global <- rep(0,k)
  cv_auc_binary <- rep(0,k)
  cv_auc_binary_global <- rep(0,k)
  warnings_auc <- rep(0,k)
  data<-data[sample(nrow(data)),]
  cv_accuracy_bin <- rep(0,k)
  cv_accuracy_balanced_bin <- rep(0,k)
  cv_f1_bin <- rep(0,k)
  cv_sens_bin <- rep(0,k)
  cv_spec_bin <- rep(0,k)


  # create fold-indices
  folds <- cut(seq(1,nrow(data)),breaks=k,labels=FALSE)
  data$fold <- folds
  for (i in 1:k){
    train <- filter(data, fold != i)
    test <- filter(data, fold == i)
    fitted_mod <- mblogit(formula, data = train, random = ~ 1|betrieb)
    pred_probs <- predict(fitted_mod, type = "response", newdata=test)
    pred_class <- factor(colnames(pred_probs)[apply(pred_probs,1,which.max)],
                         levels = c("1", "2", "3", "4", "5"))
    cm <- confusionMatrix(pred_class, test[,target])
    cv_accuarcy[i] <- cm$overall[[1]]
    cv_accuarcy_balanced[i] <- sum(cm$byClass[,6], na.rm = T)/5
    cv_macro_f1[i] <- 2 * (sum(cm$byClass[, 5], na.rm = T)/5 * sum(cm$byClass[, 6], na.rm = T)/5)/(5/sum(cm$byClass[, 5], na.rm = T) + 5/sum(cm$byClass[,6], NA.rm = T))
    cv_kappa[i] <- cm$overall[[2]]
    cv_auc[i] <- multiclass.roc(test[,target], pred_probs)$auc[1]

    ### Wsl. fuer binaere Zielgroeße antibiotikagabe hinzufuegen
    pred_probs <- as.data.frame(pred_probs)
    pred_probs$prob_antibiotikagabe <- rowSums(pred_probs[3:5])
    pred_probs$antibiotikagabe <- ifelse(pred_class %in% c(1,2), "nein", "ja")
    cv_auc_binary[i] <- auc(roc(test$antibiotikagabe, pred_probs$prob_antibiotikagabe))
    cm_binaer <- confusionMatrix(as.factor(pred_probs$antibiotikagabe), test$antibiotikagabe)
    cv_accuracy_bin[i] <- cm_binaer$overall[[1]]
    cv_accuracy_balanced_bin[i] <- cm_binaer$byClass[[11]]
    cv_f1_bin[i] <- cm_binaer$byClass[[7]]
    cv_sens_bin[i] <- cm_binaer$byClass[[1]]
    cv_spec_bin[i] <- cm_binaer$byClass[[2]]

    # Globale Schätzer
    test_global <- c(test[,target], lw_et_merged_af$verlaufsklasse)
    test_global_antibiotikagabe <- c(test$antibiotikagabe, lw_et_merged_af$antibiotikagabe)
    pred_probs_global <- add_column(pred_probs, "6" = rep(0, nrow(pred_probs)), .after = "5")
    pred_probs_global <- rbind(pred_probs_global, pred_af)
    cv_auc_global[i] <- multiclass.roc(test_global, pred_probs_global[,  c("1","2", "3", "4","5", "6")])$auc[1]
    cv_auc_binary_global[i] <- auc(roc(test_global_antibiotikagabe, pred_probs_global$prob_antibiotikagabe))
  }
  c(Reduce(paste, deparse(formula)), mean(cv_auc, na.rm = T), mean(cv_auc_global, na.rm = T),
    mean(cv_auc_binary, na.rm = T), mean(cv_auc_binary_global, na.rm = T), mean(cv_macro_f1, na.rm = T),
    mean(cv_accuarcy_balanced, na.rm = T), mean(cv_accuarcy, NA.rm = T), mean(cv_kappa, NA.rm = T),
    mean(cv_accuracy_bin, na.rm = T), mean(cv_accuracy_balanced_bin, na.rm = T),
    mean(cv_f1_bin, na.rm = T), mean(cv_sens_bin, na.rm = T), mean(cv_spec_bin, na.rm = T))
}



### Modell 1: Bestes Modell ohne Betriebseffekt
# ---------------------------- Modell und Formel laden
#source("modellwahl_multinomial")
load("overview_multi_final.RData")
final_model_multi_formula <- formula(df_overview_ordered_auc[1, "formula"])

# ---------------------------- Bestes Modell schätzen (Referenzkodierung)
final_model_multi <-
  multinom(
    final_model_multi_formula,
    data = lw_et_merged)

performance_final_model_multi <- data.frame("formula" = c(8), "cv_auc" = c(8) , "cv_auc_global" = c(8),"cv_auc_binaer" = c(8), "cv_auc_binaer_global" = c(8),
                                            "cv_macro_f1"= c(8), "cv_accuarcy_balanced"= c(8),
                                            "cv_accuarcy"= c(8), "cv_kappa"= c(8),
                                            "cv_accuracy_bin" = c(8), "cv_accuracy_balanced_bin" = c(8),
                                            "cv_f1_bin" = c(8), "cv_sens_bin" = c(8), "cv_spec_bin" = c(8))
performance_final_model_multi[1,] <- cross_validation_multinom_final(10, final_model_multi_formula, lw_et_merged, "verlaufsklasse")
performance_final_model_multi

# ---------------------------- Effectstars für bestes Modell
### Referenzkodierung
stars_ref <- star.nominal(formula = final_model_multi_formula, symmetric = FALSE, select = 2:5, nlines = 1,
                          dist.cov = 1.5, printpvalues = FALSE, test.rel = FALSE, data = lw_et_merged, refLevel = 1)
starsodds_ref <- as.data.frame(stars_ref$odds)
starsodds_ref %>%
  gt(rownames_to_stub = TRUE, caption = "Exponierte Parameterwerte") %>%
  gtsave(filename = "starodds_ref.png")

### Effektkodierung
stars_eff <- star.nominal(formula = final_model_multi_formula, symmetric = TRUE, select = 2:5, nlines = 1,
                          dist.cov = 1.5, printpvalues = FALSE, test.rel = FALSE, data = lw_et_merged)
starsodds_eff <- as.data.frame(stars_eff$odds)
starsodds_eff %>%
  gt(rownames_to_stub = TRUE, caption = "Exponierte Parameterwerte") %>%
  gtsave(filename = "starodds_eff.png")


# -----------------------------------------------------------------------------
### Modell 2: Modell mit Betriebseffekt
final_model_multi_rnd_eff <- mblogit(final_model_multi_formula, data = lw_et_merged, random = ~ 1|betrieb)
performance_final_model_multi_rnd_eff <- data.frame("formula" = c(8), "cv_auc" = c(8) , "cv_auc_global" = c(8),"cv_auc_binaer" = c(8), "cv_auc_binaer_global" = c(8),
                                                    "cv_macro_f1"= c(8), "cv_accuarcy_balanced"= c(8),
                                                    "cv_accuarcy"= c(8), "cv_kappa"= c(8),
                                                    "cv_accuracy_bin" = c(8), "cv_accuracy_balanced_bin" = c(8),
                                                    "cv_f1_bin" = c(8), "cv_sens_bin" = c(8), "cv_spec_bin" = c(8))
performance_final_model_multi_rnd_eff[1,] <- cross_validation_multinom_rnd_effect_final(10, final_model_multi_formula, lw_et_merged, "verlaufsklasse")
performance_final_model_multi_rnd_eff

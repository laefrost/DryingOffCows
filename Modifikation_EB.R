source("preprocessing_final.R")
library(mlr3)
library(mlr3learners)
library(mlr3viz)
library(dplyr)
library(rpart.plot)
library(caret)
library(pROC)
library(ROCR)
library(rpart)
library(ggparty)
set.seed(123)

### Ansatz 1: Zielgröße "Mit Auffälligkeiten vs. Unauffällig"
### -------------------------------- Kostenmatrix
costMatrix_AF <- matrix(c(0,1,10,0), ncol = 2)
rownames(costMatrix_AF) <- levels(lw_et_merged_with_all$mit_auffaelligkeiten)
colnames(costMatrix_AF) <- levels(lw_et_merged_with_all$mit_auffaelligkeiten)

### -------------------------------- Baum fitten
formula_AF <- formula(mit_auffaelligkeiten ~ laktationsnummer + zellzahl + schalmtest + schalmtest_stufe_3 + jahreszeit +
                        mean_milch + diff_milch + last_milch + trend_zz + var_milch)

tree_AF_all <- rpart(formula = formula_AF, method = "class", data = lw_et_merged_with_all, maxdepth = 4,
                     parms = list(loss=costMatrix_AF, split = "gini"))
probabilities_AF_all <- predict(tree_AF_all, type = "prob")
prediction_AF_all <- predict(tree_AF_all, type = "class")
cM_AF_all <- confusionMatrix(prediction_AF_all, lw_et_merged_with_all$mit_auffaelligkeiten)

### -------------------------------- Baum printen
rpart.plot(tree_AF_all, main ="KKB zur Identifikation von Modifikationen - Zielgröße \"Auffälligkeiten auf dem Baum\"",  roundint = FALSE, extra = 1, box.palette= list('plum3','gold2'), fallen.leaves = T)
legend(x = 0, y = 1.07, legend = c("Ja", "Nein"),
       fill = c("plum3", "gold2"),
       title = "Klassifikation in den Blättern", xpd = TRUE, bty = "o", xjust = 0)
legend(x= 0.73, y = 1.07,legend = c("Linker Wert: Anzahl \"ja\"", "Rechter Wert: Anzahl \"nein\""), #fill = NA,
       title = "Aufbau Knoten", xpd = TRUE,  bty = "o",  x.intersp = 0,  title.cex = 1, trace=FALSE)

# -----------------------------------------------------------------------------
##### Ansatz 2: Zielgröße "MBU Ergebnis minor/ohne vs. MBU Ergebnis major"
### -------------------------------- Neue ZG erstellen
lw_et_merged_with_all <- lw_et_merged_with_all %>%
  mutate(MBU_AB = case_when(MBU == "major" ~ "major",
                            TRUE ~ "minor/ohne Befund"))
lw_et_merged_with_all$MBU_AB <- as.factor(lw_et_merged_with_all$MBU_AB)
str(lw_et_merged_with_all)
### -------------------------------- Kostenmatrix
# cost sensitive decision tree, Fehler FN wiegt 10 mal so schwer
costMatrix_MBU_AB <- matrix(c(0,1,10,0), ncol = 2)
rownames(costMatrix_MBU_AB) <- levels(lw_et_merged_with_all$MBU_AB)
colnames(costMatrix_MBU_AB) <- levels(lw_et_merged_with_all$MBU_AB)

### -------------------------------- Baum fitten
formula_MBU_AB <- formula(MBU_AB ~ laktationsnummer + zellzahl + schalmtest + schalmtest_stufe_3 + jahreszeit +
                            mean_milch + diff_milch + last_milch + trend_zz + var_milch)
tree_MBU_AB <- rpart(formula = formula_MBU_AB, method = "class", data = lw_et_merged_with_all, maxdepth = 4,
                     parms = list(loss=costMatrix_AF, split = "gini"))
probabilities_MBU_AB <- predict(tree_MBU_AB, type = "prob")
prediction_MBU_AB <- predict(tree_MBU_AB, type = "class")
cM_MBU_AB <- confusionMatrix(prediction_MBU_AB, lw_et_merged_with_all$MBU_AB)

### -------------------------------- Baum printen
rpart.plot(tree_MBU_AB, main ="KKB zur Identifikation von Modifikationen - Zielgröße \"Ergebnis der MBU\"",  roundint = FALSE, extra = 1, box.palette= list('plum3','gold2'), fallen.leaves = T)
legend(x= 0.62, y=1.07, legend = c("MBU positiv", "MBU negativ/Ohne Befund"),
       fill = c("plum3", "gold2"),
       title = "Klassifikation in den Blättern", xpd = TRUE, bty = "o", xjust = 0)
legend(x= 0.62, y=0.9 ,legend = c("Linker Wert: Anzahl \"major\"", "Rechter Wert: Anzahl \"minor\\ohne Befund\""), #fill = NA,
       title = "Aufbau Knoten", xpd = TRUE,  bty = "o",  x.intersp = 0,  title.cex = 1, trace=FALSE)

# -----------------------------------------------------------------------------
### Ansatz 3: Zielgröße AB vs. kein AB
### -------------------------------- Cost matrix
costMatrix_AB <- matrix(c(0,1,100,0), ncol = 2)
rownames(costMatrix_AB) <- levels(lw_et_merged_with_all$AB)
colnames(costMatrix_AB) <- levels(lw_et_merged_with_all$AB)

### -------------------------------- Baum fitten
formula_AB <- formula(antibiotikagabe ~ laktationsnummer + zellzahl + schalmtest + schalmtest_stufe_3 + jahreszeit +
                        mean_milch + diff_milch + last_milch + trend_zz + var_milch)
tree_AB <- rpart(formula = formula_AB, method = "class", data = lw_et_merged_with_all, maxdepth = 4,
                 parms = list(loss=costMatrix_AB, split = "gini"))

### -------------------------------- Baum printen
rpart.plot(tree_AB, main ="KKB zur Identifikation von Modifikationen - Zielgröße \"Antibiotikagabe\"", roundint = FALSE, extra = 1, box.palette= list('plum3','gold2'))
legend(x= 0, y=1.07, legend = c("Ja", "Nein"),
       fill = c("plum3", "gold2"),
       title = "Klassifikation in den Blättern", xpd = TRUE, bty = "o", xjust = 0)
legend(x= 0.73, y=1.07 ,legend = c("Linker Wert: Anzahl \"ja\"", "Rechter Wert: Anzahl \"nein\""), #fill = NA,
       title = "Aufbau Knoten", xpd = TRUE,  bty = "o",  x.intersp = 0,  title.cex = 1, trace=FALSE)
levels(lw_et_merged$antibiotikagabe)
head(lw_et_merged)

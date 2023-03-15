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
  scale_y_continuous(breaks=c(0, 50, 100, 150,200, 250, 300, 350, 400, 450)) + xlab("Verlaufsklasse") +
  theme(
    axis.text = element_text(
      size = 9,
      angle = 45,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title = element_text(size = 12, face = "bold"), legend.position="bottom"
  )

# Antibiotika vs. kein Antibiotika
ggplot(data = lw_et_merged_with_all, aes(x = reorder(antibiotikagabe, antibiotikagabe, function(x) - length(x)))) +
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
tmp$mit_auffaelligkeiten <- relevel(tmp$mit_auffaelligkeiten, ref = "ja")
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
metr_vars <- lw_et_merged %>% select(zellzahl, log_zellzahl, laktationsnummer, log_laktationsnummer, mean_milch, diff_milch,
               var_milch, last_milch, trend_zz)
cor_metr <- round(cor(metr_vars),3)
melted_corr_mat <- melt(cor_metr)
# plotting the correlation heatmap
ggplot(melted_corr_mat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white") +
  coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  #label(x = "Einflussgrößen", y = "Einflussgrößen") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Korrelations-\nkoeffizient nach \nBravais-Pearson") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1)) +
  labs(y= "Einflussgrößen", x = "Einflussgrößen", title = "Heatmap für die Korrelationen zwischen den metrischen Einflussgrößen") +
  coord_fixed()

#kable_metr <- cor_metr %>%
#  kbl(escape = F,format = "latex", col.names = linebreak(c("zellzahl","log.\nZellzahl","Lakta-\ntions-\nummer", "log.\nLaktations-\nnummer",
#                              "Milch-\ndurchschnitt", "Milch-\ndifferenz",
#                              "Mich-\nvarianz","last_milch", "Zell-\nzahl-\ntrend"), align = "c")) %>%
#  kable_styling(full_width = T)
#
#writeLines(kable_metr)

cor(as.numeric(lw_et_merged$Zellzahl.max), as.numeric(lw_et_merged$diff_milch)) # 0.1772452
cor(as.numeric(lw_et_merged$Zellzahl.max), as.numeric(lw_et_merged$mean_milch)) # 0.05622195
cor(as.numeric(lw_et_merged$Zellzahl.max), as.numeric(lw_et_merged$var_milch)) # 0.1375786
cor(as.numeric(lw_et_merged$Zellzahl.max), as.numeric(lw_et_merged$last_milch)) # -0.1437737
cor(as.numeric(lw_et_merged$Zellzahl.max), as.numeric(lw_et_merged$trend_zz)) # 0.2711048

cor(as.numeric(lw_et_merged$diff_milch), as.numeric(lw_et_merged$mean_milch)) # 0.1055532
cor(as.numeric(lw_et_merged$diff_milch), as.numeric(lw_et_merged$var_milch)) # 0.8802296
cor(as.numeric(lw_et_merged$diff_milch), as.numeric(lw_et_merged$last_milch)) # -0.1399942
cor(as.numeric(lw_et_merged$diff_milch), as.numeric(lw_et_merged$trend_zz)) # 0.09037919

cor(as.numeric(lw_et_merged$mean_milch), as.numeric(lw_et_merged$trend_zz)) # 0.1390805
cor(as.numeric(lw_et_merged$mean_milch), as.numeric(lw_et_merged$var_milch)) # 0.06540648
cor(as.numeric(lw_et_merged$mean_milch), as.numeric(lw_et_merged$last_milch)) # 0.482313

cor(as.numeric(lw_et_merged$var_milch), as.numeric(lw_et_merged$trend_zz)) # 0.07769897
cor(as.numeric(lw_et_merged$var_milch), as.numeric(lw_et_merged$last_milch)) # -0.1877654

cor(as.numeric(lw_et_merged$last_milch), as.numeric(lw_et_merged$trend_zz)) # 0.006984511

cor(log(as.numeric(lw_et_merged$Zellzahl.max)), as.numeric(lw_et_merged$Zellzahl.max)) # 0.9478319 (logisch, da der Inhalt gleicht)
cor(log(as.numeric(lw_et_merged$Zellzahl.max)), as.numeric(lw_et_merged$diff_milch)) # 0.1863215
cor(log(as.numeric(lw_et_merged$Zellzahl.max)), as.numeric(lw_et_merged$mean_milch)) # 0.07100265
cor(log(as.numeric(lw_et_merged$Zellzahl.max)), as.numeric(lw_et_merged$var_milch)) # 0.1487659
cor(log(as.numeric(lw_et_merged$Zellzahl.max)), as.numeric(lw_et_merged$last_milch)) # -0.1192832
cor(log(as.numeric(lw_et_merged$Zellzahl.max)), as.numeric(lw_et_merged$trend_zz)) # 0.2505376

cor(log(as.numeric(lw_et_merged$Laktationsnummer)), log(as.numeric(lw_et_merged$Zellzahl.max)))# 0.3502199
cor(log(as.numeric(lw_et_merged$Laktationsnummer)), as.numeric(lw_et_merged$Zellzahl.max)) # 0.3270558
cor(log(as.numeric(lw_et_merged$Laktationsnummer)), as.numeric(lw_et_merged$diff_milch)) # 0.2863297
cor(log(as.numeric(lw_et_merged$Laktationsnummer)), as.numeric(lw_et_merged$mean_milch)) # 0.3862429
cor(log(as.numeric(lw_et_merged$Laktationsnummer)), as.numeric(lw_et_merged$var_milch)) # 0.2630821
cor(log(as.numeric(lw_et_merged$Laktationsnummer)), as.numeric(lw_et_merged$last_milch)) # 0.03140391
cor(log(as.numeric(lw_et_merged$Laktationsnummer)), as.numeric(lw_et_merged$trend_zz)) # 0.1712556

cor(as.numeric(lw_et_merged$Laktationsnummer), log(as.numeric(lw_et_merged$Zellzahl.max)))# 0.32052
cor(as.numeric(lw_et_merged$Laktationsnummer), as.numeric(lw_et_merged$Zellzahl.max)) # 0.3066126
cor(as.numeric(lw_et_merged$Laktationsnummer), as.numeric(lw_et_merged$diff_milch)) # 0.2481652
cor(as.numeric(lw_et_merged$Laktationsnummer), as.numeric(lw_et_merged$mean_milch)) # 0.3408543
cor(as.numeric(lw_et_merged$Laktationsnummer), as.numeric(lw_et_merged$var_milch)) # 0.2339375
cor(as.numeric(lw_et_merged$Laktationsnummer), as.numeric(lw_et_merged$last_milch)) # 0.01755883
cor(as.numeric(lw_et_merged$Laktationsnummer), as.numeric(lw_et_merged$trend_zz)) # 0.1651684
cor(as.numeric(lw_et_merged$Laktationsnummer), log(as.numeric(lw_et_merged$Laktationsnummer))) # 0.9558737 (logisch, da der Inhalt gleicht)

# Korrelationen der metrischen Variablen mit ordinalen Variablen (Kategorie) nach Spearman:
cor_spearman <- cor(metr_vars, as.numeric(lw_et_merged$kategorie),  method = "spearman", use = "complete.obs")
kable_metr <- cor_spearman %>%
    kbl(escape = F,format = "latex") %>%
    kable_styling(full_width = F, border_left = T, border_right = T)

cor(as.numeric(lw_et_merged$zellzahl), as.numeric(lw_et_merged$kategorie),  method = "spearman", use = "complete.obs") # 0.0751479

cor(as.numeric(lw_et_merged$diff_milch), as.numeric(lw_et_merged$Kategorie),  method = "spearman", use = "complete.obs") # 0.09999473

cor(as.numeric(lw_et_merged$mean_milch), as.numeric(lw_et_merged$Kategorie),  method = "spearman", use = "complete.obs") # -0.01018799

cor(as.numeric(lw_et_merged$var_milch), as.numeric(lw_et_merged$Kategorie),  method = "spearman", use = "complete.obs") # 0.1229575

cor(as.numeric(lw_et_merged$last_milch), as.numeric(lw_et_merged$Kategorie),  method = "spearman", use = "complete.obs") # -0.06196204

cor(as.numeric(lw_et_merged$trend_zz), as.numeric(lw_et_merged$Kategorie),  method = "pearson", use = "complete.obs") #  0.1404126

cor(log(as.numeric(lw_et_merged$Zellzahl.max)), as.numeric(lw_et_merged$Kategorie),  method = "pearson", use = "complete.obs") # 0.08212324

cor(log(as.numeric(lw_et_merged$Laktationsnummer)), as.numeric(lw_et_merged$Kategorie),  method = "pearson", use = "complete.obs") # -0.0004483121

cor(as.numeric(lw_et_merged$Laktationsnummer), as.numeric(lw_et_merged$Kategorie),  method = "pearson", use = "complete.obs") # 0.01836526


# Korrelationen der nominalen Variablen (Jahreszeit) mit metrischen Variablen nach Eta-Quadrat:
etaSquared(x = lm((lw_et_merged$Zellzahl.max ~ lw_et_merged$Jahreszeit))) # 0.002054021
etaSquared(x = lm((log(lw_et_merged$Zellzahl.max) ~ lw_et_merged$Jahreszeit))) # 0.004928242
etaSquared(x = lm((lw_et_merged$diff_milch ~ lw_et_merged$Jahreszeit))) # 0.004771749
etaSquared(x = lm((lw_et_merged$mean_milch ~ lw_et_merged$Jahreszeit))) # 0.008798048
etaSquared(x = lm((lw_et_merged$var_milch ~ lw_et_merged$Jahreszeit))) #  0.005429705
etaSquared(x = lm((lw_et_merged$last_milch ~ lw_et_merged$Jahreszeit))) # 0.001054507
etaSquared(x = lm((lw_et_merged$trend_zz ~ lw_et_merged$Jahreszeit))) # 0.004297034
etaSquared(x = lm((log(lw_et_merged$Laktationsnummer) ~ lw_et_merged$Jahreszeit))) # 0.006913407
etaSquared(x = lm((lw_et_merged$Laktationsnummer ~ lw_et_merged$Jahreszeit))) # 0.004420979

eta_squared(lm(lw_et_merged$Zellzahl.max ~ lw_et_merged$Jahreszeit), partial = FALSE) # 0.00205 (wie oben)
eta_squared(lm(log(lw_et_merged$Zellzahl.max) ~ lw_et_merged$Jahreszeit), partial = FALSE) # 0.00493 (wie oben)
eta_squared(lm(lw_et_merged$diff_milch ~ lw_et_merged$Jahreszeit), partial = FALSE) # 0.00477 (wie oben)
eta_squared(lm(lw_et_merged$mean_milch ~ lw_et_merged$Jahreszeit), partial = FALSE) # 0.00880 (wie oben)
eta_squared(lm(lw_et_merged$var_milch ~ lw_et_merged$Jahreszeit), partial = FALSE) # 0.00543 (wie oben)
eta_squared(lm(lw_et_merged$last_milch ~ lw_et_merged$Jahreszeit), partial = FALSE) # 0.00105 (wie oben)
eta_squared(lm(lw_et_merged$trend_zz ~ lw_et_merged$Jahreszeit), partial = FALSE) # 0.00430 (wie oben)
eta_squared(lm(log(lw_et_merged$Laktationsnummer) ~ lw_et_merged$Jahreszeit), partial = FALSE) # 0.00691 (wie oben)
eta_squared(lm(lw_et_merged$Laktationsnummer ~ lw_et_merged$Jahreszeit), partial = FALSE) # 0.00442 (wie oben)


# Korrelationen der dichotomen Variablen (Schalmtest) mit metrischen Variablen nach Punktbiseraler Korrelation:
lw_et_merged <- rename(lw_et_merged, schalmtest = SMT.Zellzahl.TS1)
cor.test(as.numeric(lw_et_merged$schalmtest), as.numeric(lw_et_merged$Zellzahl.max)) # 0.2554367
cor.test(as.numeric(lw_et_merged$schalmtest), log(as.numeric(lw_et_merged$Zellzahl.max))) # 0.242825
cor.test(as.numeric(lw_et_merged$schalmtest), as.numeric(lw_et_merged$diff_milch)) # 0.0447964
cor.test(as.numeric(lw_et_merged$schalmtest), as.numeric(lw_et_merged$mean_milch)) # -0.01957835
cor.test(as.numeric(lw_et_merged$schalmtest), as.numeric(lw_et_merged$var_milch)) # -0.009008603
cor.test(as.numeric(lw_et_merged$schalmtest), as.numeric(lw_et_merged$last_milch)) # -0.07990956
cor.test(as.numeric(lw_et_merged$schalmtest), as.numeric(lw_et_merged$trend_zz)) # 0.09795016
cor.test(as.numeric(lw_et_merged$schalmtest), log(as.numeric(lw_et_merged$Laktationsnummer))) # 0.03126443
cor.test(as.numeric(lw_et_merged$schalmtest), as.numeric(lw_et_merged$Laktationsnummer)) # 0.01721124




############# Packages installieren
#required_packages <- c("lubridate", "tidyverse")

#check_packages <- function(package){
#  if (!require(package))
#    do.call(install.packages, list(eval(package)))
#}

#lapply(required_packages, check_packages)
library(lubridate)
library(tidyverse)


# ----------------------------  Hilfsfunktionen
create_season <- function(df, date_name) {
  df$jahreszeit <- format(df[[date_name]], "%d-%m")
  # jahreszeiten definieren für ein Jahr
  WS_interval <- seq(lubridate::ymd('2016-12-21'),lubridate::ymd('2017-03-19'), by = 'days')
  FU_interval <- seq(lubridate::ymd('2017-03-20'),lubridate::ymd('2017-06-20'), by = 'days')
  SO_interval <- seq(lubridate::ymd('2017-06-21'),lubridate::ymd('2017-09-21'), by = 'days')
  HE_interval <- seq(lubridate::ymd('2017-09-22'),lubridate::ymd('2017-12-20'), by = 'days')
  # Jahr cutten
  WS_interval <- format(WS_interval, "%d-%m")
  FU_interval <- format(FU_interval, "%d-%m")
  SO_interval <- format(SO_interval, "%d-%m")
  HE_interval <- format(HE_interval, "%d-%m")

  # jahreszeitenspalte richtig befuellen
  df$jahreszeit <- as.factor(ifelse(df$jahreszeit %in% WS_interval, "Winter",
                                    ifelse(df$jahreszeit %in% FU_interval, "Fruehling",
                                           ifelse(df$jahreszeit %in% SO_interval, "Sommer", "Herbst"))))

  df
}

get_laktationsnummer <- function(row_lw){
  kuh <- row_lw["Kuh"]
  betrieb <- row_lw["Betrieb"]
  datum <- row_lw["Datum"]
  rel_entries <- et_mlp_aggregated[et_mlp_aggregated$Betrieb == betrieb & et_mlp_aggregated$Kuh == kuh, ]
  rel_entries$diff_datum <- abs(as.numeric((difftime(datum, rel_entries$letztes_Datum))))
  min_date <- min(rel_entries$diff_datum)
  rel_laktation <- rel_entries[rel_entries$diff_datum == min_date, "Laktationsnummer"]
  return(rel_laktation$Laktationsnummer[[1]])
}

tgd_raw <- readRDS("TGD_clean.Rds")
et_mlp_raw <- readRDS("MLP_clean.Rds")
lw_raw <- readRDS("LW_clean.Rds")
lw_merge_raw <- readRDS("Entscheidungsbaum(2).Rds")

# ---------------------------- Preprocessing
### TGD Datensatz
tgd_raw["Betrieb"] <- as.factor(tgd_raw$Betrieb) # 'Betrieb' als Faktor
tgd_raw["Kuh"] <- as.factor(tgd_raw$Kuh) # 'Kuh' als Faktor
#tgd_raw["Neu"] <- as.factor(tgd_raw$Neu) # 'Neu' als Faktor
tgd_raw["Befund"] <- as.factor(tgd_raw$Befund) # 'Befund' als Faktor
tgd_raw["SMT.Zellzahl"] <- as.factor(tgd_raw$SMT.Zellzahl) # 'SMT_Zellzahl' als Faktor
tgd_raw["Datum.Probenentnahme"] <- as.Date(tgd_raw$Datum.Probenentnahme, format = "%d.%m.%Y")   # 'Datum.Probenentnahme' als Datum

tgd <- tgd_raw %>% filter(Vorbericht == "7-14d vor TS") %>% select(Betrieb, Kuh, Datum.Probenentnahme, SMT.Zellzahl, Befund)

tgd_TS <- tgd_raw %>% filter(Vorbericht == "beim TS") %>% select(Betrieb, Kuh, Datum.Probenentnahme, SMT.Zellzahl)


### Einzeltier MLP Datenstaz
str(et_mlp_raw)
et_mlp <- et_mlp_raw[c("Betrieb", "Kuh", "Datum.MLP", "ST", "Laktationstage", "Laktationsnummer", "Milch.kg", "Zellzahl")]
et_mlp <- as.data.frame(et_mlp)
et_mlp["Betrieb"] <- as.factor(et_mlp$Betrieb) # 'Betrieb' als Faktor
et_mlp["Kuh"] <- as.factor(et_mlp$Kuh) # 'Kuh' als Faktor
et_mlp["Datum.MLP"] <- as.Date(et_mlp$Datum.MLP, format = "%d.%m.%Y")   # 'Datum.MLP' als Datum
et_mlp["ST"] <- as.factor(et_mlp$ST) # 'ST' als Faktor
et_mlp$Milch.kg <- gsub(",", ".", et_mlp$Milch.kg) # Kommas durch Punkte ersetzen
et_mlp["Milch.kg"] <- as.numeric(et_mlp$Milch.kg) # 'Milch.kg' numerisch machen

# Hier die Einträge vor/nach dem Trockenstellen rausfiltern
et_mlp_aggregated <- et_mlp %>% filter(ST != "S" & ST != "T" |is.na(ST)) %>%
  arrange(Betrieb, Kuh, Laktationsnummer, Datum.MLP) %>%
  group_by(Betrieb, Kuh, Laktationsnummer) %>%
  drop_na(Zellzahl) %>%
  summarise(mean_milch = mean(Milch.kg, na.rm = T), diff_milch = max(Milch.kg) - min(Milch.kg),
            letztes_Datum = last(Datum.MLP), last_milch = last(Milch.kg), var_milch = var(Milch.kg), trend_zz = lm(Zellzahl ~ Datum.MLP)$coefficients[2])

### LW Merge Datensatz
lw_merge <- lw_merge_raw[order(lw_merge_raw$Kuh), ]
v_minor <- c("Staph. chromogenes", "Staphylokokken (KNS)", "Staph. hyicus",
             "Mischinfektion_minor", "Staph. xylosus", "Staphylok.(KNS)",
             "Staph. hämolyticus", "Staph. sciuri", "Coryneforme", "Staph. epidermidis",
             "Staph. warneri", "Staph. simulans")
lw_merge <- lw_merge %>%
  mutate(MBU = case_when(Befund %in% v_minor ~ "minor",
                         Befund %in% c("kein Mastitiserreger nachweisbar", "ohne Befund") ~ "ohne Befund",
                         TRUE ~ "major"))

lw_merge <- lw_merge %>% mutate(verlaufsklasse_old_names = case_when(Stufe1 == 1 ~ "AF",
                                                           MBU == "major" ~ "AF_UF.MAJ",
                                                           MBU == "minor" & Stufe3 == 1 ~ "AF_UF.MIN.POS",
                                                           MBU == "minor" & Stufe3 == 0 ~ "AF_UF.MIN.NEG",
                                                           Stufe1 == 0 & Stufe2 == 0 & Stufe3 == 1 ~ "AF_UF.UF.POS",
                                                           Stufe1 == 0 & Stufe2 == 0 & Stufe3 == 0 ~ "UF"))
lw_merge <- lw_merge %>% mutate(verlaufsklasse = case_when(Stufe1 == 1 ~ "6",
                                                                     MBU == "major" ~ "5",
                                                                     MBU == "minor" & Stufe3 == 1 ~ "4",
                                                                     MBU == "minor" & Stufe3 == 0 ~ "2",
                                                                     Stufe1 == 0 & Stufe2 == 0 & Stufe3 == 1 ~ "3",
                                                                     Stufe1 == 0 & Stufe2 == 0 & Stufe3 == 0 ~ "1"))

lw_merge <- lw_merge %>% mutate(antibiotikagabe = case_when(verlaufsklasse_old_names %in% c("UF", "AF_UF.MIN.NEG") ~ "nein",
                                               TRUE ~ "ja"))
lw_merge <- lw_merge %>% mutate(mit_auffaelligkeiten = case_when(verlaufsklasse_old_names == "UF" ~ "nein",
                                                       TRUE ~ "ja"))

lw_merge <- lw_merge %>% select(Betrieb, Kuh, Datum, Datum.Probenentnahme.TS1, Zellzahl.max, MBU, SMT.Zellzahl.TS1, SMT.Zellzahl.TS2, Stufe1, Stufe2, Stufe3, verlaufsklasse, verlaufsklasse_old_names, antibiotikagabe, mit_auffaelligkeiten)

lw_merge["Betrieb"] <- as.factor(lw_merge$Betrieb) # 'Betrieb' als Faktor
lw_merge["Kuh"] <- as.factor(lw_merge$Kuh) # 'Kuh' als Faktor
lw_merge["verlaufsklasse_old_names"] <- as.factor(lw_merge$verlaufsklasse_old_names) # 'Kategorie' als Faktor
lw_merge["verlaufsklasse"] <- as.factor(lw_merge$verlaufsklasse) # 'Kategorie' als Faktor
lw_merge <- create_season(lw_merge, "Datum") # jahreszeit erstellen
lw_merge$Laktationsnummer <- apply(lw_merge, 1, function(x) get_laktationsnummer(x)) # Laktationsnummer erstellen
lw_merge$SMT.Zellzahl.TS1 <- as.factor(lw_merge$SMT.Zellzahl.TS1)
lw_merge$SMT.Zellzahl.TS2 <- as.factor(lw_merge$SMT.Zellzahl.TS2)


# LW Daten und aggregierte Einzeltierdaten zusammenfügen
lw_et_merged_with_all <- merge(lw_merge, et_mlp_aggregated, by = c("Betrieb", "Kuh", "Laktationsnummer"))
lw_et_merged <- lw_et_merged_with_all %>% filter(verlaufsklasse_old_names != "AF")

 lw_et_merged$Betrieb <- as.factor(parse_number(as.character(lw_et_merged$Betrieb)))
 lw_et_merged[["log_zellzahl"]] <- log(lw_et_merged$Zellzahl.max)
 lw_et_merged[["log_laktationsnummer"]] <- log(lw_et_merged$Laktationsnummer)
 lw_et_merged <- lw_et_merged %>% mutate(kategorie = case_when(
   Betrieb %in% c(2, 4, 8, 9, 13, 14, 20) ~ 2,
   TRUE ~ 1
 ))
 lw_et_merged$verlaufsklasse <- droplevels(lw_et_merged$verlaufsklasse)
 lw_et_merged$kategorie <- as.factor(lw_et_merged$kategorie)
 lw_et_merged$trend_zz <- as.numeric(lw_et_merged$trend_zz)
 lw_et_merged$antibiotikagabe <- as.factor(lw_et_merged$antibiotikagabe)
 lw_et_merged$mit_auffaelligkeiten<- as.factor(lw_et_merged$mit_auffaelligkeiten)
 lw_et_merged$verlaufsklasse <- relevel(lw_et_merged$verlaufsklasse, ref = "1")
 lw_et_merged$verlaufsklasse <- as.factor(lw_et_merged$verlaufsklasse)
 lw_et_merged$antibiotikagabe <- relevel(lw_et_merged$antibiotikagabe, ref = "nein")
 lw_et_merged$mit_auffaelligkeiten <- relevel(lw_et_merged$mit_auffaelligkeiten, ref = "nein")
 new_names <- c(schalmtest = "SMT.Zellzahl.TS1", zellzahl = "Zellzahl.max", betrieb = "Betrieb", kuh = "Kuh",
                laktationsnummer = "Laktationsnummer", schalmtest_stufe_3 = "SMT.Zellzahl.TS2")
 lw_et_merged <- rename(lw_et_merged, all_of(new_names))

 lw_et_merged_with_all$Betrieb <- as.factor(parse_number(as.character(lw_et_merged_with_all$Betrieb)))
 lw_et_merged_with_all <- lw_et_merged_with_all %>% mutate(kategorie = case_when(
   Betrieb %in% c(2, 4, 8, 9, 13, 14, 20) ~ 2,
   TRUE ~ 1
 ))
 lw_et_merged_with_all$kategorie <- as.factor(lw_et_merged_with_all$kategorie)
 lw_et_merged_with_all[["log_zellzahl"]] <- log(lw_et_merged_with_all$Zellzahl.max)
 lw_et_merged_with_all$trend_zz <- as.numeric(lw_et_merged_with_all$trend_zz)
 lw_et_merged_with_all[["log_laktationsnummer"]] <- log(lw_et_merged_with_all$Laktationsnummer)
 lw_et_merged_with_all$verlaufsklasse <- relevel(lw_et_merged_with_all$verlaufsklasse, ref = "1")
 lw_et_merged_with_all$antibiotikagabe <- as.factor(lw_et_merged_with_all$antibiotikagabe)
 lw_et_merged_with_all$mit_auffaelligkeiten<- as.factor(lw_et_merged_with_all$mit_auffaelligkeiten)
 lw_et_merged_with_all <- rename(lw_et_merged_with_all, all_of(new_names))
 lw_et_merged_with_all$antibiotikagabe <- relevel(lw_et_merged_with_all$antibiotikagabe, ref = "nein")
 lw_et_merged_with_all$mit_auffaelligkeiten <- relevel(lw_et_merged_with_all$mit_auffaelligkeiten, ref = "nein")


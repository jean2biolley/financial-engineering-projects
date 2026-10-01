
# ==============================================================================
# Code associated with the LLSMS2224 Forecasting Group Project 2025-2026
#
# Professor:  Bertrand Candelon
#
# Jean de Biolley - 55612100
# Evrard de Marchant – 51102000
# ==============================================================================

library(readr)
library(zoo)
library(ggplot2)
library(tseries)
library(mFilter)
library(forecast)
library(vars)
library(gridExtra)
library(quantmod)
library(selectr)
library(tidyverse)
library(readxl)
library(dplyr)
library(scales)
library(tidyr)
library(lmtest)
library(dlm)

library(imputeTS)
library(urca)
library(forecast)

###############################################################################
###                                                                         ###
###     Question 1 : Economic theory of the 10 years bond yield spread      ###              
###                                                                         ###
###############################################################################
### Import from FRED 

# Définition des symboles et des noms
symbols <- c("IRLTLT01ESQ156N", "IRLTLT01DEQ156N", "CP0000ESM086NEST", 
             "CP0000DEM086NEST", "LRHUTTTTESQ156S", "LRHUTTTTDEQ156S", 
             "CPMNACSCAB1GQES", "CPMNACSCAB1GQDE","XTNTVA01ESQ664S", "XTNTVA01DEQ664S") 

names_map <- c("GovBYSP", "GovBYGE", "HICPSP", "HICPGE", "URSP", "URGE", 
               "gdpSP", "gdpGE","TBCommoditiesSP", "TBCommoditiesGE")

# Création d'un environnement pour stocker les données
env <- new.env()

# Importation des données
getSymbols(
  Symbols = symbols,
  src = "FRED",
  from = "2013-10-01",
  to   = "2025-05-31",
  env  = env
)

# Fusion des données importées
raw <- do.call(merge, eapply(env, function(x) x))

# Mapping correct FRED -> Nouveaux noms
mapping <- c(
  "IRLTLT01ESQ156N" = "GovBYSP", 
  "IRLTLT01DEQ156N" = "GovBYGE", 
  "CP0000ESM086NEST" = "HICPSP",
  "CP0000DEM086NEST" = "HICPGE",
  "LRHUTTTTESQ156S" = "URSP",
  "LRHUTTTTDEQ156S" = "URGE",
  "CPMNACSCAB1GQES" = "gdpSP",
  "CPMNACSCAB1GQDE" = "gdpGE",
  "XTNTVA01ESQ664S" = "TBCommoditiesSP",
  "XTNTVA01DEQ664S" = "TBCommoditiesGE"
)

# Renommer les colonnes selon le mapping
colnames(raw) <- mapping[colnames(raw)]

####### --- Transformation en Data Frame ---

df <- data.frame(
  Date = as.yearqtr(index(raw)),
  GovBYSP = as.numeric(raw$GovBYSP),
  GovBYGE = as.numeric(raw$GovBYGE),
  HICPSP = as.numeric(raw$HICPSP),
  HICPGE = as.numeric(raw$HICPGE),
  URSP = as.numeric(raw$URSP),
  URGE = as.numeric(raw$URGE),
  gdpSP = as.numeric(raw$gdpSP),
  gdpGE = as.numeric(raw$gdpGE),
  TBCommoditiesSP = as.numeric(raw$TBCommoditiesSP),
  TBCommoditiesGE = as.numeric(raw$TBCommoditiesGE)
)

# Aggregate monthly to quarterly mean (CPI & TS)
df_q2013 <- df |>
  group_by(Date) |>
  summarize(
    GovBYSP  = mean(GovBYSP,  na.rm = TRUE),
    GovBYGE  = mean(GovBYGE,  na.rm = TRUE),
    HICPSP   = mean(HICPSP,   na.rm = TRUE),
    HICPGE   = mean(HICPGE,   na.rm = TRUE),
    URSP     = mean(URSP,     na.rm = TRUE),
    URGE     = mean(URGE,     na.rm = TRUE),
    gdpSP    = mean(gdpSP,    na.rm = TRUE),
    gdpGE    = mean(gdpGE,    na.rm = TRUE),
    TBCommoditiesSP = mean(TBCommoditiesSP,    na.rm = TRUE),
    TBCommoditiesGE = mean(TBCommoditiesGE,    na.rm = TRUE)
  )

df_q2013 <- na.omit(df_q2013)
#jai a partir de 01/10/2013 et dcp je calcule vite l'Inf et puis jenelve cette date la pour navoir que des 
#q1 2014
df_q2013$InfSP <- c(NA, diff(log(df_q2013$HICPSP)))
df_q2013$InfGE <- c(NA, diff(log(df_q2013$HICPGE)))
df_q <- na.omit(df_q2013)
# Inflation (Log-différence du HICP) car je ne peux pas plus tard faire le difflog pour les growth rates 
# sur le spread de HICP SP et GE car je calculerai le taux de croissance de l'écart
# et ce n'est pas une mesure standard. Ce qu'on veut en économie, c'est le Différentiel d'Inflation (Inflation Espagne - Inflation Allemagne).
#
df_q$TBCommoditiesSP <- df_q$TBCommoditiesSP / 1000000000
df_q$TBCommoditiesGE <- df_q$TBCommoditiesGE / 1000000000
#
#
df_q$gdpSP <- df_q$gdpSP / 1000
df_q$gdpGE <- df_q$gdpGE / 1000
#
#
#transfo de GovBYSP ET GE en bps et de InfSP ET GE en %
#et pas en bps car en part 3 je dois lui faire 2 diff et il est enorme sinon par rapport aux autres 
#(sinon jaurai du lui faire 2 fois *100 car de base ils
# sont en decimales 0,02 (in one hundreth of percent) ->puis je fais *100 ca fait 2 pour lavoir en %
#-> puis je fais encore *100 ca fait 200 pour lavoir en bps car gemini ma dit de faire comme ca
df_q$GovBYSP <- df_q$GovBYSP * 100
df_q$GovBYGE <- df_q$GovBYGE * 100
df_q$InfSP <- df_q$InfSP * 100
df_q$InfGE <- df_q$InfGE * 100
#
######################### vision rapides des variables sans leurs spread #####################
#
p_GovBY <- ggplot(df_q, aes(x = Date)) +
  geom_line(aes(y = GovBYSP, color = "Spain")) +
  geom_line(aes(y = GovBYGE, color = "Germany")) +
  labs(title = "Gov Bond Yield", y = "Gov Bond Yield (in bps)", x = "Year") +
  theme_minimal()
print(p_GovBY)

#
p_Inf <- ggplot(df_q, aes(x = Date)) +
  geom_line(aes(y = InfSP, color = "Spain")) +
  geom_line(aes(y = InfGE, color = "Germany")) +
  labs(title = "Inflation", y = "Inflation (in %)", x = "Year") +
  theme_minimal()
print(p_Inf)

#
p_UR <- ggplot(df_q, aes(x = Date)) +
  geom_line(aes(y = URSP, color = "Spain")) +
  geom_line(aes(y = URGE, color = "Germany")) +
  labs(title = "Unemployment Rate (%)", y = "Unemployment Rate (%)", x = "Year") +
  theme_minimal()
print(p_UR)

#
p_gdp <- ggplot(df_q, aes(x = Date)) +
  # Ajout de la ligne à 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_line(aes(y = gdpSP, color = "Spain")) +
  geom_line(aes(y = gdpGE, color = "Germany")) +
  labs(title = "Gross Domestic Product", y = "Gross Domestic Product (in B of Euros)", x = "Year") +
  theme_minimal()
print(p_gdp)

#
p_TBCommodities <- ggplot(df_q, aes(x = Date)) +
  geom_line(aes(y = TBCommoditiesSP, color = "Spain")) +
  geom_line(aes(y = TBCommoditiesGE, color = "Germany")) +
  labs(title = "Trade Balance Commodities", y = "Trade Balance Commodities (in B of Euros)", x = "Year") +
  theme_minimal()
print(p_TBCommodities)

############## fin vision rapide des variables sans leurs spread ###########

###############import de 4 time series ds lexcel:  ##################

#CDS
CDS <- read_excel("variables pour forecast.xlsx", 
                  sheet = "Combined CDS", col_types = c("date", 
                                                        "numeric", "numeric"))

# Renommer proprement les colonnes
CDS <- CDS |>
  rename(
    CDSSpain   = `Spain CDS; Quarterly`,
    CDSGermany = `Germany CDS; Quarterly`
  )
# 1. Transformer la date de l'Excel en format "yearqtr" (trimestre)
CDS$Date <- as.yearqtr(CDS$Date)
# On repart de ton fichier importé et on cree un CDS_clean qui est exactement comme df_q en forme dcp facilement left_joignable
CDS_clean <- CDS |>
  mutate(Date = as.yearqtr(Date)) |> # Convertit les dates mensuelles en trimestres
  group_by(Date) |>
  summarize(
    CDSSpain = mean(CDSSpain, na.rm = TRUE),
    CDSGermany = mean(CDSGermany, na.rm = TRUE)
  )
# On fusionne sur le tableau trimestriel propre
df_q <- df_q |>
  left_join(CDS_clean, by = "Date")




#Debt to gdp

Debt <- read_excel("variables pour forecast.xlsx", 
                   sheet = "debt to gdp ", col_types = c("date", 
                                                         "numeric", "numeric"))
# Renommer proprement les colonnes
Debt <- Debt |>
  rename(
    DebtSpain   = `Government debt Spain (consolidated), percent of GDP`,
    DebtGermany = `Government debt Germany (consolidated), percent of GDP`
  )
# 1. Transformer la date de l'Excel en format "yearqtr" (trimestre) cad "2013 Q4"
Debt$Date <- as.yearqtr(Debt$DATE)
# On repart de ton fichier importé et on cree un CDS_clean qui est exactement comme df_q en forme dcp facilement left_joignable
Debt_clean <- Debt |>
  mutate(Date = as.yearqtr(Date)) |> # Convertit les dates mensuelles en trimestres
  group_by(Date) |>
  summarize(
    DebtSpain = mean(DebtSpain, na.rm = TRUE),
    DebtGermany = mean(DebtGermany, na.rm = TRUE)
  )
# On fusionne sur le tableau trimestriel propre
df_q <- df_q |>
  left_join(Debt_clean, by = "Date")


#Budget Balance per GDP
Budgetbalance <- read_excel("variables pour forecast.xlsx", 
                            sheet = " Budget Balance (% GDP)", col_types = c("date", 
                                                                             "numeric", "numeric"))
Budgetbalance <- Budgetbalance |>
  rename(
    BudgetbalanceSpain   = `Spain budget Balance; USD; Q`,
    BudgetbalanceGermany = `Germany Budgert Balance; USD; Q`,
  )
# 1. Transformer la date de l'Excel en format "yearqtr" (trimestre) cad "2013 Q4"
Budgetbalance$Date <- as.yearqtr(Budgetbalance$Date)
# On repart de ton fichier importé et on cree un CDS_clean qui est exactement comme df_q en forme dcp facilement left_joignable
Budgetbalance_clean <- Budgetbalance |>
  mutate(Date = as.yearqtr(Date)) |> # Convertit les dates mensuelles en trimestres
  group_by(Date) |>
  summarize(
    BudgetbalanceSpain = mean(BudgetbalanceSpain, na.rm = TRUE),
    BudgetbalanceGermany = mean(BudgetbalanceGermany, na.rm = TRUE)
  )
# On fusionne sur le tableau trimestriel propre
df_q <- df_q |>
  left_join(Budgetbalance_clean, by = "Date")





#VSTOXX et il est mtn en quaterly apres le run de ce petit code ci dessous
Vstoxx <- read_excel("variables pour forecast.xlsx", 
                     sheet = "STOXX 50 Volatility VSTOXX EUR ", 
                     col_types = c("date", "numeric"))
# Renommer proprement les colonnes
Vstoxx <- Vstoxx |>
  rename(
    VstoxxEUROPE   = `Dernier`,
  )
# --- CORRECTION DE LA VIRGULE ---
# On remplace "," par "." et on convertit en numérique
# gsub(ce_qu'on_cherche, ce_par_quoi_on_remplace, la_colonne)
Vstoxx$VstoxxEUROPE <- as.numeric(gsub(",", ".", Vstoxx$VstoxxEUROPE))
# 1. Transformer la date de l'Excel en format "yearqtr" (trimestre) cad "2013 Q4"
Vstoxx$Date <- as.yearqtr(Vstoxx$Date)
# On repart de ton fichier importé et on cree un CDS_clean qui est exactement comme df_q en forme dcp facilement left_joignable
Vstoxx_clean <- Vstoxx |>
  mutate(Date = as.yearqtr(Date)) |> # Convertit les dates mensuelles en trimestres
  group_by(Date) |>
  summarize(
    VstoxxEUROPE = mean(VstoxxEUROPE, na.rm = TRUE)
  )
# On fusionne sur le tableau trimestriel propre
df_q <- df_q |>
  left_join(Vstoxx_clean, by = "Date")
##################### fin import ds lexcel###############


######################### vision rapides des variables excels sans leurs spread #####################
#
p_CDS <- ggplot(df_q, aes(x = Date)) +
  geom_line(aes(y = CDSSpain, color = "Spain")) +
  geom_line(aes(y = CDSGermany, color = "Germany")) +
  labs(title = "CDS (bps)", y = "CDS (bps)", x = "Year") +
  theme_minimal()
print(p_CDS)

#
p_Debt <- ggplot(df_q, aes(x = Date)) +
  geom_line(aes(y = DebtSpain, color = "Spain")) +
  geom_line(aes(y = DebtGermany, color = "Germany")) +
  labs(title = "Debt per gdp (%)", y = "Debt per gdp (%)", x = "Year") +
  theme_minimal()
print(p_Debt)

#
p_Budgetbalance <- ggplot(df_q, aes(x = Date)) +
  geom_line(aes(y = BudgetbalanceSpain, color = "Spain")) +
  geom_line(aes(y = BudgetbalanceGermany, color = "Germany")) +
  labs(title = "Budget balance (%)", y = "Budget balance (%)", x = "Year") +
  theme_minimal()
print(p_Budgetbalance)

#
#
p_Vstoxx <- ggplot(df_q, aes(x = Date)) +
  geom_line(aes(y = VstoxxEUROPE, color = "Europe")) +
  labs(title = "VSTOXX (volatility of the EUROSTOXX50 Index)", y = "(in %)", x = "Year") +
  theme_minimal()
print(p_Vstoxx)
############## fin vision rapide des variables excels sans leurs spread ###########


############################ Calcul de tout les spreads ######################

# 3. Calcul des Spreads et que faire Spain - Germany meme si ca va faire des nombres negatifs m'a dit gemini
# Spread GovBY (Spain - Germany)
df_q$SpreadGovBY <- (df_q$GovBYSP - df_q$GovBYGE)
# Spread Inf 
df_q$SpreadInf <- (df_q$InfSP - df_q$InfGE)
# Spread UR
df_q$SpreadUR  <- (df_q$URSP - df_q$URGE)
# Spread GDP 
df_q$SpreadGDP <- (df_q$gdpSP - df_q$gdpGE)
# Spread TBCommodities 
df_q$SpreadTBCommodities <- (df_q$TBCommoditiesSP - df_q$TBCommoditiesGE)
# Spread CDS
df_q$SpreadCDS <- (df_q$CDSSpain - df_q$CDSGermany)
# Spread Debt per gdp
df_q$SpreadDebt <- (df_q$DebtSpain - df_q$DebtGermany)
# Spread Budgetbalance per gdp
df_q$SpreadBudgetbalance <- (df_q$BudgetbalanceSpain - df_q$BudgetbalanceGermany)
#pas de spread pour Vstoxx dcp car il represente la volatility des 50 actions Europeennes
############################ fin Calcul de tout les spreads ######################

############## vision rapide de SpreadGovBY SP-GE en partliculier pour le report ###########
p_SpreadGovBY <- ggplot(df_q, aes(x = Date, y = SpreadGovBY)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "Spread Gov Bond Yield Spain-Germany (2014-2025)", x = "Date", y = "(in bps)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
grid.arrange(p_SpreadGovBY, ncol = 1)
############## fin vision rapide de SpreadGovBY SP-GE en partliculier pour le report ###########

##################### plots de tous les spreads sur un seul plot ############
# P1: Spread GovBY
p1 <- ggplot(df_q, aes(x = Date, y = SpreadGovBY)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "Spread Gov Bond Yield (2014-2025)", x = "Date", y = "(in bps)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

# P2: Spread Inf
p2 <- ggplot(df_q, aes(x = Date, y = SpreadInf)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(color = "red", linewidth = 1) +
  labs(title = "Spread Inf (2014-2025)", x = "Date", y = "(in %)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

# P3: Spread UR
p3 <- ggplot(df_q, aes(x = Date, y = SpreadUR)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(color = "darkgreen", linewidth = 1) +
  labs(title = "Spread Unemployment Rate (2014-2025)", x = "Date", y = "(in %)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

# P4: Spread GDP
p4 <- ggplot(df_q, aes(x = Date, y = SpreadGDP)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(color = "orange", linewidth = 1) +
  labs(title = "Spread GDP (2014-2025)", x = "Date", y = "(in B of Euros)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

# P5: Spread TBCommodities
p5 <- ggplot(df_q, aes(x = Date, y = SpreadTBCommodities)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(color = "gold", linewidth = 1) +
  labs(title = "Spread TradeBalance Commodities (2014-2025)", x = "Date", y = "(in B of Euros)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

# P6: Spread CDS
p6 <- ggplot(df_q, aes(x = Date, y = SpreadCDS)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(color = "lightblue", linewidth = 1) +
  labs(title = "Spread CDS (2014-2025)", x = "Date", y = "(in bps)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

# P7: Spread Debt
p7 <- ggplot(df_q, aes(x = Date, y = SpreadDebt)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(color = "pink", linewidth = 1) +
  labs(title = "Spread Debt / GDP (2014-2025)", x = "Date", y = "(in % of GDP)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

# P8: Spread Budget Balance 
p8 <- ggplot(df_q, aes(x = Date, y = SpreadBudgetbalance)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(color = "purple", linewidth = 1) +
  labs(title = "Spread Budgetbalance / GDP (2014-2025)", x = "Date", y = "(in % of GDP)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

# P9: VSTOXX Europe
p9 <- ggplot(df_q, aes(x = Date, y = VstoxxEUROPE)) +
  # Pas de ligne à 0 car la volatilité est toujours > 0
  geom_line(color = "brown", linewidth = 1) +
  labs(title = "Vstoxx Europe (2014-2025)", x = "Date", y = "(in %)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

# Affichage final avec les 9 graphiques
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)

################### les histograms de mes variables #########
# 1. Création du dataframe "df_histo" (Format Large)
df_histo <- data.frame(
  Date = df_q$Date,
  SpreadGovBY = df_q$SpreadGovBY,
  SpreadInf = df_q$SpreadInf,
  SpreadUR = df_q$SpreadUR,
  SpreadGDP = df_q$SpreadGDP,
  SpreadTBCommodities = df_q$SpreadTBCommodities,
  SpreadCDS = df_q$SpreadCDS,
  SpreadDebt = df_q$SpreadDebt,
  SpreadBudgetbalance = df_q$SpreadBudgetbalance,
  VstoxxEUROPE = df_q$VstoxxEUROPE
)
# 2. Transformation en format "Long" cad tout en 1 colonne pour les 9 variables toutes les 2014 Q1 lun a la suite de lautre (pour les graphiques groupés)
df_long <- df_histo |>
  pivot_longer(
    cols = c(SpreadGovBY, SpreadInf, SpreadUR, SpreadGDP, 
             SpreadTBCommodities, SpreadCDS, SpreadDebt, SpreadBudgetbalance,VstoxxEUROPE),
    names_to = "Variable",
    values_to = "Value"
  )
#plot des histogrammes
# On définit l'option globale pour aider  (sécurité anti-scientifique)
options(scipen = 999)
p_hist <- ggplot(df_long, aes(x = Value, fill = Variable)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~ Variable, scales = "free") +
  labs(
    title = "Histograms of Variables",
    x = "Values",
    y = "Frequency",
    fill = "Variable"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  # C'est ici que j'ai changé le code pour corriger l'erreur :
  scale_x_continuous(labels = function(x) format(x, big.mark = " ", scientific = FALSE))
print(p_hist)
###################### fin les histograms de mes variables #########

########### faire les First Differences et les plots ###############

# 2. Calcul des First Differences et non des Growth Rates et 
#La première différence te montre si l'écart s'est creusé entre SP-GE (valeur positive) ou réduit entre SP-GE (valeur négative) 
#par rapport au trimestre précédent (=l’observation precedente).
#->mtn df_histo a tt les variables et leurs firstdiff en pas cleaned en outliers !!! (voir prt 3 et part 4 debut)
df_histo$diffSpreadGovBY <- c(NA, diff(df_histo$SpreadGovBY))
df_histo$diffSpreadInf <- c(NA, diff(df_histo$SpreadInf))
df_histo$diffSpreadUR <- c(NA, diff(df_histo$SpreadUR))
df_histo$diffSpreadGDP <- c(NA, diff(df_histo$SpreadGDP))
df_histo$diffSpreadTBCommodities <- c(NA, diff(df_histo$SpreadTBCommodities))
df_histo$diffSpreadCDS <- c(NA, diff(df_histo$SpreadCDS))
df_histo$diffSpreadDebt <- c(NA, diff(df_histo$SpreadDebt))
df_histo$diffSpreadBudgetbalance <- c(NA, diff(df_histo$SpreadBudgetbalance))
df_histo$diffVstoxxEUROPE <- c(NA, diff(df_histo$VstoxxEUROPE))
# 4. Format Long pour le graphique
df_long <- df_histo %>%
  pivot_longer(cols = starts_with("diff"), 
               names_to = "Variable", 
               values_to = "Value")
# 6. Plot combiné (Tout sur le même graphique)
ggplot(df_long, aes(x = Date, y = Value, color = Variable)) +
  geom_line(linewidth = 1) +
  labs(title = "First Differences of All Variables (2014-2025)",
       x = "Date",
       y = "Change",
       color = "Variables" # Titre de la légende
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    # --- CHANGEMENT 2 : ON AFFICHE LA LÉGENDE ---
    legend.position = "right", # Met la légende à droite (ou "bottom" pour en bas)
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_yearqtr(format = "%Y Q%q", n = 10) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c(
    "diffSpreadGovBY" = "blue",
    "diffSpreadInf" = "red",
    "diffSpreadUR" = "darkgreen",
    "diffSpreadGDP" = "orange",
    "diffSpreadTBCommodities" = "lightblue",
    "diffSpreadCDS" = "gold",
    "diffSpreadDebt" = "pink",
    "diffSpreadBudgetbalance" = "purple",
    "diffVstoxxEUROPE" = "brown"
  ))

##### Histograms of First Differences 
# On définit l'option globale pour aider (sécurité anti-scientifique)
options(scipen = 999)
# On utilise le df_long qui contient maintenant les First Differences
p_hist_diff <- ggplot(df_long, aes(x = Value, fill = Variable)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  # scales = "free" est toujours crucial ici car les unités sont différentes
  facet_wrap(~ Variable, scales = "free") + 
  labs(
    title = "Histograms of First Differences of Variables", # Titre mis à jour
    x = "Values (Change)",
    y = "Frequency",
    fill = "Variable"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 8), # Taille réduite car tes noms de variables sont longs
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12),
    legend.position = "none" # Pas besoin de légende grâce aux titres des facettes
  ) +
  # Ton formatage anti-scientifique
  scale_x_continuous(labels = function(x) format(x, big.mark = " ", scientific = FALSE))
print(p_hist_diff)
###### fin Histograms of First Differences 
############ fin First Differences et leurs plots ##############

###############################################################################
###                                                                         ###
###                 Question 2 : Univariate model ARIMA                     ###
###                                                                         ###
###############################################################################
#recreer une nouvelle petite df jiuste pour spreadGovBY
#gemini me dit de faire le exploratory data analysis (EDA) sur le spread en niveau cad bps et puis 
#sur le first difference du spread (cad que larima aura deja d=1)
#EDA sur le spread en level:
# 1. Création du Dataframe avec la valeur BRUTE (Niveau)
df_SpreadGovBY_Level <- data.frame(
  Date = df_q$Date, 
  SpreadGovBY_Level = df_q$SpreadGovBY # Pas de diff, pas de log, juste la valeur brute
)
# On nettoie les NA éventuels (issus du merge par exemple)
df_SpreadGovBY_Level <- na.omit(df_SpreadGovBY_Level)
# --- P10: Line Plot (Time Series) ---
p10 <- ggplot(df_SpreadGovBY_Level, aes(x = Date, y = SpreadGovBY_Level)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  # Titre adapté : On précise que c'est le "Level"
  labs(title = "Spread Gov Bond Yield - Level (2014-2025)", 
       x = "Date", 
       y = "Spread (in bps)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_yearqtr(format = "%Y Q%q", n = 10) +
  scale_y_continuous(labels = comma)
# --- P11: Histogramme ---
p11 <- ggplot(df_SpreadGovBY_Level, aes(x = SpreadGovBY_Level)) +
  geom_histogram(bins = 30, fill = "darkgreen", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Spread Gov Bond Yield",
       x = "Spread (in bps)", 
       y = "Frequency") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(labels = comma)
# --- P12: Boxplot ---
p12 <- ggplot(df_SpreadGovBY_Level, aes(y = SpreadGovBY_Level)) +
  geom_boxplot(fill = "darkgreen", outlier.color = "red", outlier.shape = 16, alpha = 0.7) +
  labs(title = "Boxplot of Spread Gov Bond Yield",
       y = "Spread (in bps)", 
       x = "") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank()) +
  scale_y_continuous(labels = comma)
# --- Affichage combiné ---
grid.arrange(p10, p11, p12, ncol = 1)


# --- P13: Line Plot (Après traitement des outliers) ---
p13 <- ggplot(df_SpreadGovBY_Level, aes(x = Date, y = SpreadGovBY_Level)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  labs(title = "Spread Gov Bond Yield - Level (Outliers Treated)", 
       x = "Date", 
       y = "Spread (in bps)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_yearqtr(format = "%Y Q%q", n = 10) +
  scale_y_continuous(labels = comma)
# --- P14: Histogramme (Après traitement des outliers) ---
p14 <- ggplot(df_SpreadGovBY_Level, aes(x = SpreadGovBY_Level)) +
  geom_histogram(bins = 30, fill = "darkgreen", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Spread (Outliers Treated)",
       x = "Spread (in bps)", 
       y = "Frequency") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(labels = comma)
# --- P15: Boxplot (Après traitement des outliers) ---
p15 <- ggplot(df_SpreadGovBY_Level, aes(y = SpreadGovBY_Level)) +
  geom_boxplot(fill = "darkgreen", outlier.color = "red", outlier.shape = 16, alpha = 0.7) +
  labs(title = "Boxplot of Spread (Outliers Treated)",
       y = "Spread (in bps)", 
       x = "") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank()) +
  scale_y_continuous(labels = comma)

grid.arrange(p13, p14, p15, ncol = 1)

# --- Test de Stationnarité (ADF) ---
# On teste la stationnarité sur la série nettoyée des outliers
adf_result <- adf.test(df_SpreadGovBY_Level$SpreadGovBY_Level)
print(adf_result)
# Interprétation automatique pour t'aider
if(adf_result$p.value < 0.05) {
  cat("P-value < 0.05 : We reject H0. The series is STATIONARY.\n")
} else {
  cat("P-value >= 0.05 : We fail to reject H0. The series is NON-STATIONARY (Differencing required).\n")
}
#fin EDA sur le spread en level

#EDA sur le spread en first diff:
#le test adf dit quil faut encore diff dcp faire la first difference de SpreadGivBY
# 1. Création du Dataframe spécifique pour l'analyse des différences
df_SpreadGovBY_Diff <- data.frame(
  Date = df_q$Date, 
  SpreadGovBY_Diff = df_q$SpreadGovBY 
)
# 2. Création de la variable "firstdiff"
#en fait cisdessus la colonne "SpreadGovBY_Diff" nest pas encore en first diff car on lui a juste mit "df_q$SpreadGovBY " dedans
# mais apres on fait first diff calcul juste cidessous:
df_SpreadGovBY_Diff$SpreadGovBY_Diff <- c(NA, diff(df_SpreadGovBY_Diff$SpreadGovBY_Diff))
# 3. On nettoie les NA (la première ligne sera supprimée à cause du diff) et dcp on etait encore en q1 2014 et 
#on sera plus que en q2 2014
df_SpreadGovBY_Diff <- na.omit(df_SpreadGovBY_Diff)
# --- P16: Line Plot (First Difference) ---
p16 <- ggplot(df_SpreadGovBY_Diff, aes(x = Date, y = SpreadGovBY_Diff)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "Spread Gov Bond Yield - First Difference", 
       x = "Date", 
       y = "Change (in bps)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_yearqtr(format = "%Y Q%q", n = 10) +
  scale_y_continuous(labels = comma)
# --- P17: Histogramme (First Difference) ---
p17 <- ggplot(df_SpreadGovBY_Diff, aes(x = SpreadGovBY_Diff)) +
  geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of First Difference",
       x = "Change (in bps)", 
       y = "Frequency") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(labels = comma)
# --- P18: Boxplot (First Difference) ---
p18 <- ggplot(df_SpreadGovBY_Diff, aes(y = SpreadGovBY_Diff)) +
  geom_boxplot(fill = "blue", outlier.color = "red", outlier.shape = 16, alpha = 0.7) +
  labs(title = "Boxplot of First Difference",
       y = "Change (in bps)", 
       x = "") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank()) +
  scale_y_continuous(labels = comma)
grid.arrange(p16, p17, p18, ncol = 1)

Q1_diff <- quantile(df_SpreadGovBY_Diff$SpreadGovBY_Diff, 0.25)  
Q3_diff <- quantile(df_SpreadGovBY_Diff$SpreadGovBY_Diff, 0.75)
IQR_diff <- Q3_diff - Q1_diff  
lower_bound_diff <- Q1_diff - 1.5 * IQR_diff
upper_bound_diff <- Q3_diff + 1.5 * IQR_diff

outliers_diff <- df_SpreadGovBY_Diff %>%
  filter(SpreadGovBY_Diff < lower_bound_diff | SpreadGovBY_Diff > upper_bound_diff)
non_outlier_mean_diff <- df_SpreadGovBY_Diff %>%
  filter(SpreadGovBY_Diff >= lower_bound_diff & SpreadGovBY_Diff <= upper_bound_diff) %>%
  summarise(mean_value = mean(SpreadGovBY_Diff)) %>%
  pull(mean_value)
# Vision
print(outliers_diff)
print(non_outlier_mean_diff)
#
df_SpreadGovBY_Diff <- df_SpreadGovBY_Diff %>%
  mutate(SpreadGovBY_Diff = ifelse(SpreadGovBY_Diff < lower_bound_diff | SpreadGovBY_Diff > upper_bound_diff, 
                                   non_outlier_mean_diff, SpreadGovBY_Diff))

# --- P19: Line Plot (Après traitement) ---
p19 <- ggplot(df_SpreadGovBY_Diff, aes(x = Date, y = SpreadGovBY_Diff)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "Spread First Difference (Outliers Treated)", 
       x = "Date", 
       y = "Change (in bps)") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_yearqtr(format = "%Y Q%q", n = 10) +
  scale_y_continuous(labels = comma)
# --- P20: Histogramme (Après traitement) ---
p20 <- ggplot(df_SpreadGovBY_Diff, aes(x = SpreadGovBY_Diff)) +
  geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of First Difference (Outliers Treated)",
       x = "Change (in bps)", 
       y = "Frequency") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(labels = comma)
# --- P21: Boxplot (Après traitement) ---
p21 <- ggplot(df_SpreadGovBY_Diff, aes(y = SpreadGovBY_Diff)) +
  geom_boxplot(fill = "blue", outlier.color = "red", outlier.shape = 16, alpha = 0.7) +
  labs(title = "Boxplot of First Difference (Outliers Treated)",
       y = "Change (in bps)", 
       x = "") +
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank()) +
  scale_y_continuous(labels = comma)
grid.arrange(p19, p20, p21, ncol = 1)

# On teste la stationnarité sur la série différenciée et nettoyée
adf_result_diff <- adf.test(df_SpreadGovBY_Diff$SpreadGovBY_Diff)
print(adf_result_diff)
# Interprétation automatique
if(adf_result_diff$p.value < 0.05) {
  cat("P-value < 0.05 : We reject H0. The series is STATIONARY.\n")
} else {
  cat("P-value >= 0.05 : We fail to reject H0. The series is NON-STATIONARY (Differencing required).\n")
}
# fin EDA sur le spread en first diff

#debut Arima
# 1. ACF & PACF Plots
# Permet de visualiser les autocorrélations pour estimer les ordres p et q
par(mfrow = c(1, 2))
Acf(df_SpreadGovBY_Diff$SpreadGovBY_Diff, main = "ACF of Spread First Diff")
Pacf(df_SpreadGovBY_Diff$SpreadGovBY_Diff, main = "PACF of Spread First Diff")
# Reset de l'affichage graphique
par(mfrow = c(1, 1)) 
# stationary = TRUE car on utilise déjà la série différenciée et les 3 autres criteres seasonal,
#,... font juste que le aut.arima ets plus precis
# 2. Estimation du modèle ARIMA
# Comme je vois pas clairement sur ACF et PACF que les barres tombes apres x lags, le mieux cest dutiliser 
#aut.arima qui va choisir pour moi
#jai fait larima sur q2 2014 - q2 2025 
arima_model <- auto.arima(df_SpreadGovBY_Diff$SpreadGovBY_Diff, 
                          seasonal = FALSE, 
                          stationary = TRUE, 
                          stepwise = FALSE, 
                          include.drift= TRUE, 
                          approximation = FALSE)
summary(arima_model) #il me chosiit ARIMA (0,0,2) cad p=0 et q=2 dcp MA(2) qui veut dire que lobs actuelle de mon spread 
#en first difference est influencé par les erreurs des 2 trimestres precedents a chaque fois

# 3. Récupération des valeurs ajustées (Fitted Values)
df_SpreadGovBY_Diff$fit <- fitted(arima_model)
# 4. Plot Actual vs Fitted
ggplot(df_SpreadGovBY_Diff, aes(x = Date)) +
  # Actual values (Bleu continu)
  geom_line(aes(y = SpreadGovBY_Diff, color = "Actual Values"), linewidth = 1) +  
  # Fitted values (Rouge pointillé)
  geom_line(aes(y = fit, color = "Fitted Values"), linewidth = 1, linetype = "dashed") +  
  labs(
    title = "Actual vs. Fitted Values of ARIMA Model",
    x = "Date",
    y = "Change (in bps)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  scale_x_yearqtr(format = "%Y Q%q", n = 10) +
  scale_color_manual(values = c("Actual Values" = "blue", "Fitted Values" = "red"))
#fin ARIMA

###        Question 2 (Suite) : Diagnostics & Forecast (2026-2028)          ###
# 5. Extraction et Analyse des Résidus Globaux
residualsarima <- residuals(arima_model)
# Test de Normalité (Shapiro-Wilk)
# H0: Les résidus suivent une loi normale. Si p-value > 0.05, on ne rejette pas H0 (C'est bon).
shapiro_test <- shapiro.test(residualsarima) 
# Test d'Hétéroscédasticité (Breusch-Pagan)
# H0: Variance constante (Homoscédasticité). Si p-value > 0.05, c'est bon.
# On teste les résidus par rapport aux valeurs ajustées (fitted values)
bp_test <- bptest(residualsarima ~ df_SpreadGovBY_Diff$fit) #le ~ veut dire je pense "on va chercehr les residuals de fit ds df_SpreadGovBY_Diff
# --- AJOUT : Test d'Autocorrélation (Breusch-Godfrey) ---
# H0: Pas d'autocorrélation jusqu'à l'ordre spécifié. Si p-value > 0.05, c'est bon (Bruit blanc).
# order = 4 car données trimestrielles (on vérifie la mémoire sur 1 an/4 trimestres).
# "~ 1" signifie qu'on teste les résidus contre une constante (intercept) et leurs propres lags.
bg_test_global <- bgtest(residualsarima ~ 1, order = 4)
# Calcul des erreurs globales
mae <- mean(abs(residualsarima)) #comme un test sur qualite des residuals de larima
rmse <- sqrt(mean(residualsarima^2)) #comme un test sur qualite des residuals de larima

# Affichage des résultats globaux
cat("\n--- Global Model Diagnostics ---\n")
print(shapiro_test)
print(bp_test)
print(bg_test_global) # Affichage du Breusch-Godfrey
cat("Global MAE:", mae, "\n")
cat("Global RMSE:", rmse, "\n")


### --- FORECAST 2026-2028 (Avec intervalles 80% et 95%) --- ###
# 1. Génération des dates futures (2026-2028)
last_date <- max(df_SpreadGovBY_Diff$Date)
future_dates <- last_date + seq(0.25, by = 0.25, length.out = 12)
# 2. Calcul du forecast (sur 12 périodes)
forecast_result <- forecast(arima_model, h = 12)
# 3. Création du dataframe avec les DEUX intervalles
# [, 1] correspond au 80% et [, 2] au 95%
df_forecastarima_plot <- data.frame(
  Date = future_dates,
  Value = as.numeric(forecast_result$mean),
  Lower80 = as.numeric(forecast_result$lower[, 1]), 
  Upper80 = as.numeric(forecast_result$upper[, 1]),
  Lower95 = as.numeric(forecast_result$lower[, 2]), 
  Upper95 = as.numeric(forecast_result$upper[, 2]),
  Type = "Forecast"
)
# 4. Plot manuel avec ggplot2
plot_forecast <- ggplot() +
  # --- COUCHE 1 : Données historiques ---
  geom_line(data = df_SpreadGovBY_Diff, 
            aes(x = Date, y = SpreadGovBY_Diff, color = "Historical Data"), 
            linewidth = 1) +
  # --- COUCHE 2 : Intervalle 95% (Zone plus claire en dessous) ---
  geom_ribbon(data = df_forecastarima_plot, 
              aes(x = Date, ymin = Lower95, ymax = Upper95, fill = "95% Confidence Interval"), 
              alpha = 0.2) +
  # --- COUCHE 3 : Intervalle 80% (Zone plus foncée au-dessus) ---
  geom_ribbon(data = df_forecastarima_plot, 
              aes(x = Date, ymin = Lower80, ymax = Upper80, fill = "80% Confidence Interval"), 
              alpha = 0.3) +
  # --- COUCHE 4 : Ligne du Forecast ---
  geom_line(data = df_forecastarima_plot, 
            aes(x = Date, y = Value, color = "Forecast"), 
            linewidth = 1, linetype = "dashed") +
  # --- ESTHÉTIQUE & LABELS ---
  labs(
    title = "Forecast ARIMA: Spread Gov Bond Yield (First Difference) 2026-2028",
    x = "Date",
    y = "Change of the Spread (in bps)",
    color = "Legend",
    fill = "Confidence Intervals"
  ) +
  theme_gray() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  # --- GESTION DES COULEURS ---
  scale_color_manual(values = c("Historical Data" = "black", "Forecast" = "blue")) +
  # On utilise deux nuances de bleu pour différencier les deux zones
  scale_fill_manual(values = c("95% Confidence Interval" = "skyblue3", 
                               "80% Confidence Interval" = "steelblue")) +
  # --- GESTION DE L'AXE X ---
  scale_x_yearqtr(format = "%Y Q%q", n = 15)
print(plot_forecast)


### --- Evaluation Pre-Covid vs Post-Covid (Corrigé avec Breusch-Godfrey) --- ###
# 1. On ajoute les résidus au dataframe principal
# Note : Assurez-vous que 'arima_model' est bien votre objet modèle créé juste avant
df_SpreadGovBY_Diff$res <- residuals(arima_model) 
# 2. Création des sous-ensembles (Split)
# Pre-Covid : Avant 2020
pre_covid <- subset(df_SpreadGovBY_Diff, Date < 2020)
# Post-Covid : A partir de 2020 inclus
post_covid <- subset(df_SpreadGovBY_Diff, Date >= 2020)
# 3. Calcul des métriques d'erreur par période
mae_pre <- mean(abs(pre_covid$res))
rmse_pre <- sqrt(mean(pre_covid$res^2))
mae_post <- mean(abs(post_covid$res))
rmse_post <- sqrt(mean(post_covid$res^2))
# 4. Tests Statistiques par période
# --- Période PRE-COVID ---
shapiro_test_pre <- shapiro.test(pre_covid$res) 
# Breusch-Godfrey Test
# H0 : Pas d'autocorrélation en série jusqu'à l'ordre spécifié.
# order = 4 car données trimestrielles (on vérifie l'année complète) et modèle MA(2) inclus dedans.
bg_test_pre <- bgtest(pre_covid$res ~ 1, order = 4) 
# Breusch-Pagan Test (Hétéroscédasticité)
bp_pre <- bptest(pre_covid$res ~ pre_covid$fit)
# --- Période POST-COVID ---
shapiro_test_post <- shapiro.test(post_covid$res) 
# Breusch-Godfrey Test (order = 4)
bg_test_post <- bgtest(post_covid$res ~ 1, order = 4)
# Breusch-Pagan Test
bp_post <- bptest(post_covid$res ~ post_covid$fit)
# 5. Affichage comparatif des résultats
cat("\n###############################################################\n")
cat("###        COMPARISON PRE-COVID vs POST-COVID RESULTS       ###\n")
cat("###############################################################\n\n")
cat("--- ERREURS (Accuracy) ---\n")
cat(sprintf("MAE  Pre-Covid: %.4f  |  MAE  Post-Covid: %.4f\n", mae_pre, mae_post))
cat(sprintf("RMSE Pre-Covid: %.4f  |  RMSE Post-Covid: %.4f\n", rmse_pre, rmse_post))
cat("\n--- NORMALITÉ (Shapiro-Wilk) ---\n")
cat("Pre-Covid p-value: ", shapiro_test_pre$p.value, "\n")
cat("Post-Covid p-value:", shapiro_test_post$p.value, "\n")
cat("\n--- AUTOCORRÉLATION (Breusch-Godfrey, lag=4) ---\n")
cat("Pre-Covid p-value: ", bg_test_pre$p.value, "\n")
cat("Post-Covid p-value:", bg_test_post$p.value, "\n")
cat("\n--- HÉTÉROSCÉDASTICITÉ (Breusch-Pagan) ---\n")
cat("Pre-Covid p-value: ", bp_pre$p.value, "\n")
cat("Post-Covid p-value:", bp_post$p.value, "\n")



############### Backtesting 2020-2022
# 1. Séparation des données (Splitting)
# Train : On s'arrête juste avant 2020 (fin 2019)
df_train_bt <- subset(df_SpreadGovBY_Diff, Date < 2020) #subset refais le meme df genre meme meme noms de colonnes et tt mais juste prend une partie des dates
#je lappelle "train" ici car pour le backtesting, le modele doit sentrainer sur 2014-2019
# Test (Réalité) : On prend les 3 années suivantes (12 trimestres) pour comparer
# Note : On prend exactement 12 trimestres à partir de 2020 pour matcher le h=12
df_test_bt <- subset(df_SpreadGovBY_Diff, Date >= 2020)#ce subset refais aussi le meme df genre meme meme noms de colonnes et tt mais juste prend une partie des dates
#je lappelle "test" ici car pour le backtesting, le modele doit faire un test de forecast cad forecaster pour 2020-2022
df_test_bt <- head(df_test_bt, 12) # On garde juste 2020, 2021, 2022
# 2. Estimation du modèle UNIQUEMENT sur le passé (2014-2019)
arima_backtest <- auto.arima(df_train_bt$SpreadGovBY_Diff, 
                             seasonal = FALSE, 
                             stationary = TRUE, 
                             stepwise = FALSE, 
                             approximation = FALSE)
cat("\n--- Backtest Model Summary (Trained on 2014-2019) ---\n")
summary(arima_backtest) #ici on trouve un arima (0,0,0) mais cest normal cest que la variable firsft diff spread
#avec son observation actuelle netait pas influencée par des erruers de trimestres precedents
#et je ne peux pas utiliser mon arima (0,0,2) car il a ete fait sur aussi 2020-2025 et je commets ce 
#qu'on appelle un biais de survie ou de la "triche" statistique : ton modèle utilise des paramètres 
#influencés par une crise qu'il n'est pas censé connaître encore.
# 3. Forecast sur 3 ans (12 trimestres) "Out-of-Sample"
forecast_bt <- forecast(arima_backtest, h = 12)
# 4. Création du Dataframe manuel pour le plot (Comme demandé)
# Génération des dates futures pour le backtest (2020 Q1 -> 2022 Q4)
last_date_train <- max(df_train_bt$Date)
dates_bt <- last_date_train + seq(0.25, by = 0.25, length.out = 12)
df_forecast_bt_plot <- data.frame(
  Date = dates_bt,
  Value = as.numeric(forecast_bt$mean),
  # Extraction des 80% (colonne 1)
  Lower80 = as.numeric(forecast_bt$lower[, 1]), 
  Upper80 = as.numeric(forecast_bt$upper[, 1]),
  # Extraction des 95% (colonne 2)
  Lower95 = as.numeric(forecast_bt$lower[, 2]), 
  Upper95 = as.numeric(forecast_bt$upper[, 2]),
  Type = "Forecast Backtest"
)
# 5. Calcul des erreurs du Backtest (Forecast vs Réalité)
# On compare la colonne 'Value' du forecast avec la colonne 'SpreadGovBY_Diff' de df_test_bt
error_bt <- df_test_bt$SpreadGovBY_Diff - df_forecast_bt_plot$Value #les vraies obs - les valeurs du 
#forecast (arima(0,0,0)) dcp ca fait que des 0
mae_bt_score <- mean(abs(error_bt))
rmse_bt_score <- sqrt(mean(error_bt^2))
cat("\n###############################################################\n")
cat("###                 BACKTESTING RESULTS                     ###\n")
cat("###############################################################\n")
cat(sprintf("Backtest MAE  (2020-2022): %.4f\n", mae_bt_score))
cat(sprintf("Backtest RMSE (2020-2022): %.4f\n", rmse_bt_score))
# 6. Plot du Backtest : Histoire (Noir) + Réalité (Rouge) + Forecast (Bleu)
plot_backtest <- ggplot() +
  # --- COUCHE 1 : Les données d'entraînement (2014-2019) ---
  geom_line(data = df_train_bt, 
            aes(x = Date, y = SpreadGovBY_Diff, color = "Training Data (Pre-2020)"), 
            linewidth = 1) +
  # --- COUCHE 2 : La Réalité Observée (2020-2022) ---
  geom_line(data = df_test_bt, 
            aes(x = Date, y = SpreadGovBY_Diff, color = "Actual Observed (2020-2022)"), 
            linewidth = 1) +
  # --- COUCHE 3 : Intervalle 95% (Zone externe - plus claire) ---
  geom_ribbon(data = df_forecast_bt_plot, 
              aes(x = Date, ymin = Lower95, ymax = Upper95, fill = "95% Confidence Interval"), 
              alpha = 0.3) +
  # --- COUCHE 4 : Intervalle 80% (Zone interne - plus foncée) ---
  geom_ribbon(data = df_forecast_bt_plot, 
              aes(x = Date, ymin = Lower80, ymax = Upper80, fill = "80% Confidence Interval"), 
              alpha = 0.5) +
  # --- COUCHE 5 : La ligne du Forecast Backtest ---
  geom_line(data = df_forecast_bt_plot, 
            aes(x = Date, y = Value, color = "Forecast (Backtest)"), 
            linewidth = 1, linetype = "dashed") +
  # --- ESTHÉTIQUE ---
  labs(
    title = "Backtesting ARIMA: Forecast vs Actuals (Covid Period 2020-2022)",
    subtitle = paste("Model trained on 2014-2019 data only. RMSE:", round(rmse_bt_score, 2)),
    x = "Date",
    y = "Change of the Spread (in bps)",
    color = "Legend",
    fill = "Confidence Interval"
  ) +
  theme_gray() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  # --- COULEURS ---
  scale_color_manual(values = c(
    "Training Data (Pre-2020)" = "black", 
    "Actual Observed (2020-2022)" = "red", 
    "Forecast (Backtest)" = "blue"
  )) +
  # Définition des couleurs demandées pour les intervalles
  scale_fill_manual(values = c("95% Confidence Interval" = "skyblue3", 
                               "80% Confidence Interval" = "steelblue")) +
  # --- AXE X ---
  scale_x_yearqtr(format = "%Y Q%q", n = 15)
# Affichage
print(plot_backtest)
########## fin backtesting





###############################################################################
##########################       Part 3     ###############################


###############################################################################

# 1. Création du dataframe
df_histo <- na.omit(df_histo)
#df_histo reste le meme sans lui avoir remplacé ses outliers pcq cest sur df_part3 quon remplace les outliers
#et dcp df_histo pour etre reutilisé pour faire le DIAGNOSTIC : VÉRIFICATION DES OUTLIERS REMPLACÉS aussi en part 4
df_part3 <- df_histo[, c(
  "Date",
  "diffSpreadGovBY",
  "diffSpreadInf",
  "diffSpreadUR",
  "diffSpreadGDP",
  "diffSpreadTBCommodities",
  "SpreadCDS",
  "diffSpreadDebt",
  "diffSpreadBudgetbalance",
  "diffVstoxxEUROPE"
)]

######################## verifying stationarity #####################
replace_outliers_iqr <- function(data, column_name) {
  if (!column_name %in% colnames(data)) {
    stop("Specified column does not exist in the dataset")
  }
  Q1 <- quantile(data[[column_name]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column_name]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  non_outlier_mean <- data %>%
    filter(data[[column_name]] >= lower_bound & data[[column_name]] <= upper_bound) %>%
    summarise(mean_value = mean(data[[column_name]], na.rm = TRUE)) %>%
    pull(mean_value)
  data <- data %>%
    mutate(!!column_name := ifelse(
      data[[column_name]] < lower_bound | data[[column_name]] > upper_bound,
      non_outlier_mean,
      data[[column_name]]
    ))
  return(data)
}
df_part3 <- replace_outliers_iqr(df_part3, "diffSpreadGovBY")
df_part3 <- replace_outliers_iqr(df_part3, "diffSpreadInf")
df_part3 <- replace_outliers_iqr(df_part3, "diffSpreadUR")
df_part3 <- replace_outliers_iqr(df_part3, "diffSpreadGDP")
df_part3 <- replace_outliers_iqr(df_part3, "diffSpreadTBCommodities")
df_part3 <- replace_outliers_iqr(df_part3, "SpreadCDS")
df_part3 <- replace_outliers_iqr(df_part3, "diffSpreadDebt")
df_part3 <- replace_outliers_iqr(df_part3, "diffSpreadBudgetbalance")
df_part3 <- replace_outliers_iqr(df_part3, "diffVstoxxEUROPE")

# --- DIAGNOSTIC : VÉRIFICATION DES OUTLIERS REMPLACÉS ---
# 1. Liste des variables à vérifier (toutes sauf la Date)
vars_to_check <- c("diffSpreadGovBY", "diffSpreadInf", "diffSpreadUR", 
                   "diffSpreadGDP", "diffSpreadTBCommodities", "SpreadCDS", 
                   "diffSpreadDebt", "diffSpreadBudgetbalance", "diffVstoxxEUROPE")
# 2. Boucle pour comparer df_histo (Avant sans outliers enlevés) et df_part3 (Après)
cat("\n=== RAPPORT DE REMPLACEMENT DES OUTLIERS ===\n")
total_outliers <- 0

for (var in vars_to_check) {
  # On repère les indices où la valeur est différente
  # (On utilise round pour éviter les minuscules erreurs de virgule flottante)
  changed_idx <- which(round(df_histo[[var]], 6) != round(df_part3[[var]], 6))
  
  if (length(changed_idx) > 0) {
    cat(paste0("\nVariable : ", var, " -> ", length(changed_idx), " modifications.\n"))
    
    # On affiche le détail pour cette variable
    for (i in changed_idx) {
      date_val <- df_part3$Date[i]
      old_val <- df_histo[[var]][i]
      new_val <- df_part3[[var]][i]
      
      cat(sprintf("  - %s : Avant = %.4f  ->  Après = %.4f\n", 
                  as.character(date_val), old_val, new_val))
    }
    total_outliers <- total_outliers + length(changed_idx)
  } else {
    cat(paste0("\nVariable : ", var, " -> Aucun outlier détecté.\n"))
  }
}
cat(paste0("\nTotal d'outliers remplacés dans le dataset : ", total_outliers, "\n"))
cat("============================================\n")

adf_results <- sapply(df_part3 %>% dplyr::select(-Date), function(x) {
  tryCatch({
    test <- adf.test(na.omit(x))
    return(test$p.value)
  }, error = function(e) return(NA))
})
adf_summary <- data.frame(
  Series = names(adf_results),
  P_Value = adf_results,
  Stationary = ifelse(adf_results < 0.05, "Yes", "No")
)
print(adf_summary)

#gemini me dit que pour VARX, je dois rediff une 2e fois les non statio et garder en first statio les autres
#mais pour le vecemx on peut prendre des non statio
#dcp ceux qui ne sont pas statio en firstduff cest: diffSpreadGDP, diffSpreadTBCommodities, diffSpreadDebt,
#diffSpreadBudgetbalance, diffVstoxxEUROPE

# 2. On applique la 2e différence UNIQUEMENT sur les variables non-stationnaires
df_part3$diff2SpreadGDP           <- c(NA, diff(df_part3$diffSpreadGDP))
df_part3$diff2SpreadTBCommodities <- c(NA, diff(df_part3$diffSpreadTBCommodities))
df_part3$diff2SpreadDebt          <- c(NA, diff(df_part3$diffSpreadDebt))
df_part3$diff2SpreadBudgetbalance        <- c(NA, diff(df_part3$diffSpreadBudgetbalance))
df_part3$diff2VstoxxEUROPE              <- c(NA, diff(df_part3$diffVstoxxEUROPE))

adf_results <- sapply(df_part3 %>% dplyr::select(-Date), function(x) {
  tryCatch({
    test <- adf.test(na.omit(x))
    return(test$p.value)
  }, error = function(e) return(NA))
})
adf_summary <- data.frame(
  Series = names(adf_results),
  P_Value = adf_results,
  Stationary = ifelse(adf_results < 0.05, "Yes", "No")
)
print(adf_summary)
# 3. On nettoie la ligne de NA créée par le nouveau diff()
df_part3 <- na.omit(df_part3)

# 4. On fait un tri des colonnes ds df_part3 et on sélectionne les bonnes colonnes pour le VARX
# On prend les "vertes" (1ère diff) et les nouvelles "rouges" (2ème diff)
df_part3 <- df_part3[, c("Date", 
                         "diffSpreadGovBY", "diffSpreadInf", "diffSpreadUR","diff2SpreadGDP",
                         "diff2SpreadTBCommodities", "diffSpreadCDS", "diff2SpreadDebt", 
                         "diff2SpreadBudgetbalance", "diff2VstoxxEUROPE")]

# 1. Préparation des données pour le VARselect
# On retire la colonne "Date" car VARselect ne travaille que sur des données numériques
df_varx_for_selection <- df_part3 %>% 
  dplyr::select(-Date)

# 2. Exécution de VARselect
# Nous testons jusqu'à 10 lags 
# 'type = const' inclut une constante dans le modèle
lag_selection <- VARselect(df_varx_for_selection, lag.max = 10, type = "const")
# 3. Affichage des résultats
cat("\n=== CRITÈRES DE SÉLECTION DU NOMBRE DE LAGS ===\n")
print(lag_selection$selection)
#il a chosi 4 lags max

#
granger_test_bidirectional <- function(data, var1, var2, max_lag) {
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data frame.")
  }
  # We use lapply to loop through each lag from 1 to max_lag
  results <- lapply(1:max_lag, function(lag) {
    
    # Direction 1: Does var2 cause var1?
    test1 <- grangertest(as.formula(paste(var1, "~", var2)), order = lag, data = data)
    
    # Direction 2: Does var1 cause var2?
    test2 <- grangertest(as.formula(paste(var2, "~", var1)), order = lag, data = data)
    
    # We combine both results into a small data frame
    data.frame(
      Lag_Order = rep(lag, 2),
      Test_Direction = c(paste(var2, "causes", var1), paste(var1, "causes", var2)),
      F_Statistic = c(test1$F[2], test2$F[2]),
      P_Value = c(test1$`Pr(>F)`[2], test2$`Pr(>F)`[2])
    )
  })
  # We bind all the rows together to get one clear table
  results_df <- do.call(rbind, results)
  return(results_df)
}
#
# on met max_lag = 4 comme on avait trouvé ci dessus comme gemini ma dit de faire
granger_test_bidirectional(df_varx_for_selection, "diffSpreadGovBY", "diffSpreadInf", 4)
granger_test_bidirectional(df_varx_for_selection, "diffSpreadGovBY", "diffSpreadUR", 4)
granger_test_bidirectional(df_varx_for_selection, "diffSpreadGovBY", "diff2SpreadGDP", 4)
granger_test_bidirectional(df_varx_for_selection, "diffSpreadGovBY", "diff2SpreadTBCommodities", 4)
granger_test_bidirectional(df_varx_for_selection, "diffSpreadGovBY", "diffSpreadCDS", 4)
granger_test_bidirectional(df_varx_for_selection, "diffSpreadGovBY", "diff2SpreadDebt", 4)
granger_test_bidirectional(df_varx_for_selection, "diffSpreadGovBY", "diff2SpreadBudgetbalance", 4)
granger_test_bidirectional(df_varx_for_selection, "diffSpreadGovBY", "diff2VstoxxEUROPE", 4)
#gemini me conseille: Regarde plutot les p-value du Lag 1 ou 2 sur les 4 au total: Si elle est significative cad <0.05 
#cad on rejette H0: "la var influence lautre var". Si on a dcp 2 rejets de H0 dans les 2 sens pour les 2 variables,
#et bien dcp on sait que la variable (pas principale) en question a une influence sur tt le systeme et est en plus 
# directe et rapide et quelle est endogene dcp. Dcp au vu des resltats il n'y a aucune variable pas principale qui
# est vrmt endogene mais en meme temps gemini me conseille de prendre CDS quand meme comme var endogene pour pouvoir 
#avoir un vrai VARX, VECMX et Kalman et pas un bete ARIMA et surtout pcq la theorie dit que nrmt CDS influence tt le
#systeme. CDS a sur ses 3 lags <0.05 cad quil influence Spread mais que dans lautre sens il na pas moins que 0.31
#comme p-value dcp le Spread ninfluence pas vrmt CDS. (On a aussi vu que le VARX et VECEMX avec que CDS comme endo
#faisait les meilleurs resultats MAE/RSME/les 3 tests que quand on prenanit UR avec ou VstoxxEUROPE avec)

# 1. Séparation des variables ENDOGÈNES et EXOGÈNES
df_var_endo <- df_part3[, c("diffSpreadGovBY", "diffSpreadCDS")]
df_var_exo <- df_part3[, c("diff2SpreadGDP", "diffSpreadUR")]
# 2. Sélection du Lag Optimal
# On informe VARselect qu'il y a des variables exogènes pour ajuster les critères d'information
lag_selectionfinalvarx <- VARselect(df_var_endo, lag.max = 4, type = "const", exogen = df_var_exo)
optimal_lagfinalvarx <- lag_selectionfinalvarx$selection["AIC(n)"]
#->le systeme a choisi AIC(n)= 3 dcp optimal_lagfinalvarx = 3
cat("\n--- CRITÈRES DE SÉLECTION DU LAG (VARX) ---\n")
print(lag_selectionfinalvarx$selection)
cat(paste("Optimal Lag chosen (AIC):", optimal_lagfinalvarx, "\n"))
# Sécurité : Si AIC suggère 0 (rare), on force à 1
if (optimal_lagfinalvarx == 0) optimal_lagfinalvarx <- 1

# 3. Estimation du modèle VARX
varx_model <- VAR(y = df_var_endo, p = optimal_lagfinalvarx, type = "const", exogen = df_var_exo)
summary(varx_model)

# 4. Extraction des valeurs ajustées et des résidus
# fitted() et residuals() retournent une matrice avec les 2 variables endogènes
fitvarx <- as.data.frame(fitted(varx_model))
residualsvarx <- as.data.frame(residuals(varx_model))

# 7. Graphical Comparison: Actual vs Fitted (GGPLOT2 Version)
# A. Préparation du Dataframe aligné
# Le modèle fitted a moins d'observations que df_part3 à cause des lags (optimal_lagfinalvarx)
# On crée un dataframe temporaire qui ne contient que les dates communes
n_fitted <- nrow(fitvarx)
df_varx_plot <- data.frame(
  Date = tail(df_part3$Date, n_fitted),           # Les dates correspondantes
  Actual = tail(df_part3$diffSpreadGovBY, n_fitted), # Les vraies valeurs correspondantes
  Fitted = fitvarx$diffSpreadGovBY                # Les valeurs estimées
)
# B. Plot Actual vs Fitted (Style identique à ARIMA)
p_varx_fit <- ggplot(df_varx_plot, aes(x = Date)) +
  # Actual values (Bleu continu)
  geom_line(aes(y = Actual, color = "Actual Values"), linewidth = 1) +  
  # Fitted values (Rouge pointillé)
  geom_line(aes(y = Fitted, color = "Fitted Values"), linewidth = 1, linetype = "dashed") +  
  labs(
    title = "VARX: Actual vs. Fitted Values (Gov Bond Yield Spread)",
    x = "Date",
    y = "Spread Variation (in bps)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  scale_x_yearqtr(format = "%Y Q%q", n = 10) +
  scale_color_manual(values = c("Actual Values" = "blue", "Fitted Values" = "red"))
print(p_varx_fit)

############forecasting de VARX
# 1. Préparation des variables exogènes pour le futur (2026-2028)
n_ahead <- 12
# Hypothèse : On utilise la dernière valeur connue des exogènes (Last Value Carried Forward)
last_exog <- tail(df_var_exo, 1)
future_exog <- last_exog[rep(1, n_ahead), ]
rownames(future_exog) <- NULL 
# 2. Exécution des Forecasts (Pour 95% ET 80%)
# On lance la prédiction deux fois pour récupérer les bornes des deux intervalles
varx_forecast_95 <- predict(varx_model, n.ahead = n_ahead, dumvar = future_exog, ci = 0.95)
varx_forecast_80 <- predict(varx_model, n.ahead = n_ahead, dumvar = future_exog, ci = 0.80)
varx_forecast_80
# 3. Création du Dataframe pour le plot
# Récupération de la dernière date historique
last_date_hist <- max(df_part3$Date)
# Création de la séquence de dates futures
future_dates <- last_date_hist + seq(0.25, by = 0.25, length.out = n_ahead)
# Extraction des données
# On cible la variable 'diffSpreadGovBY' dans la liste 'fcst'
fcst_matrix_95 <- varx_forecast_95$fcst$diffSpreadGovBY
fcst_matrix_80 <- varx_forecast_80$fcst$diffSpreadGovBY
df_forecast_varx <- data.frame(
  Date = future_dates,
  Forecast = fcst_matrix_95[, "fcst"], # La moyenne est la même pour les deux
  # Bornes 95%
  Lower95 = fcst_matrix_95[, "lower"],
  Upper95 = fcst_matrix_95[, "upper"],
  # Bornes 80%
  Lower80 = fcst_matrix_80[, "lower"],
  Upper80 = fcst_matrix_80[, "upper"]
)
# 4. Plot du Forecast (Style "Backtest" avec doubles intervalles)
p_forecast_varx <- ggplot() +
  # --- COUCHE 1 : Historique (Noir) ---
  geom_line(data = df_part3, 
            aes(x = Date, y = diffSpreadGovBY, color = "Historical"), 
            linewidth = 1) +
  # --- COUCHE 2 : Intervalle 95% (Zone externe - plus claire) ---
  geom_ribbon(data = df_forecast_varx, 
              aes(x = Date, ymin = Lower95, ymax = Upper95, fill = "95% Confidence Interval"), 
              alpha = 0.3) +
  # --- COUCHE 3 : Intervalle 80% (Zone interne - plus foncée) ---
  geom_ribbon(data = df_forecast_varx, 
              aes(x = Date, ymin = Lower80, ymax = Upper80, fill = "80% Confidence Interval"), 
              alpha = 0.5) +
  # --- COUCHE 4 : Forecast (Bleu pointillé) ---
  geom_line(data = df_forecast_varx, 
            aes(x = Date, y = Forecast, color = "Forecast"), 
            linewidth = 1, linetype = "dashed") +
  # --- ESTHÉTIQUE ---
  labs(
    title = "VARX Forecast: Spread Gov Bond Yield (2026-2028)",
    subtitle = "Exogenous variables held constant at last observed value",
    y = "Spread Variation (in bps)", 
    x = "Date",
    color = "Legend",
    fill = "Confidence Intervals"
  ) +
  theme_gray() + # Utilisation de theme_gray() comme dans ton backtest
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  # --- COULEURS ---
  scale_color_manual(values = c(
    "Historical" = "black", 
    "Forecast" = "blue"
  )) +
  # --- COULEURS DE REMPLISSAGE (FILL) ---
  scale_fill_manual(values = c(
    "95% Confidence Interval" = "skyblue3", 
    "80% Confidence Interval" = "steelblue"
  )) + 
  # --- FORMAT DE L'AXE X ---
  scale_x_yearqtr(format = "%Y Q%q", n = 10)
print(p_forecast_varx)

# 5. Impulse Response Function (IRF) : utilise Bootstrap avec 100 runs
# Analyse de choc : Comment le Spread réagit à un choc sur les CDS ?
#engros ce plot cest avec les modele mathematique que jai calculé ci dessus, sur la time series du CDS je simule un 
#choc genre d'un coup il va tres fort en haut et le plot bootstrap me montre comme le spreadGovBY va rajuster le tir
#genre au debut le Spread va trop haut ds lestimation du CDS apres 2 periodes (cad 2x3 mois cad 6 mois)
#puis il va trop bas ds lesti du CDS mais deja un peu moins fort apres 4 periodes cad apres 12 mois et apres 6 periodes cad 
#1,5 ans bah il arrive presque totalement a etre proche du CDS et a bien le predire
#cest pour prouver que en effet le CDS influence bien le Spread (par contre comme on a vu, linverse est faux selon
# la granger causality results)
irf_cds <- irf(varx_model, impulse = "diffSpreadCDS", response = "diffSpreadGovBY", n.ahead = 8, boot = TRUE)
# Affichage des IRF
plot(irf_cds, main = "Shock: CDS -> Spread")
par(mfrow = c(1, 1)) # Remettre la fenêtre graphique normale


####################Diagnostics du VARX
shapiro_testvarx <- shapiro.test(residualsvarx$diffSpreadGovBY)
bg_testvarx <- bgtest(residualsvarx$diffSpreadGovBY ~ 1, order = 4)
bp_testvarx <- bptest(residualsvarx$diffSpreadGovBY ~ fitvarx$diffSpreadGovBY)
maevarx <- mean(abs(residualsvarx$diffSpreadGovBY))
rmsevarx <- sqrt(mean(residualsvarx$diffSpreadGovBY^2))


# 6. Affichage des résultats Diagnostics
cat("\n--- DIAGNOSTICS POUR LA CIBLE (diffSpreadGovBY) ---\n")
cat(sprintf("Shapiro-Wilk Normality Test (p-value): %.4f\n", shapiro_testvarx$p.value))
cat(sprintf("Breusch-Godfrey Serial Correlation Test (p-value): %.4f\n", bg_testvarx$p.value))
cat(sprintf("Breusch-Pagan Heteroskedasticity Test (p-value): %.4f\n", bp_testvarx$p.value))
cat("--------------------------------------------------\n")
cat(sprintf("Mean Absolute Error (MAE): %.6f\n", maevarx))
cat(sprintf("Root Mean Squared Error (RMSE): %.6f\n", rmsevarx))

############# 6. Evaluation Pre/Post Covid (Comme demandé)
# On crée un DF temporaire avec résidus et fit pour l'analyse
df_varxeval <- data.frame(
  Date = df_part3$Date[(optimal_lagfinalvarx + 1):nrow(df_part3)], # Ajustement des dates perdues par les lags
  Res = residualsvarx$diffSpreadGovBY,
  Fit = fitvarx$diffSpreadGovBY
)
# Split Pre/Post Covid
pre_covid_varx <- subset(df_varxeval, Date < 2020)
post_covid_varx <- subset(df_varxeval, Date >= 2020)
# Calcul des métriques
mae_pre <- mean(abs(pre_covid_varx$Res))
rmse_pre <- sqrt(mean(pre_covid_varx$Res^2))
mae_post <- mean(abs(post_covid_varx$Res))
rmse_post <- sqrt(mean(post_covid_varx$Res^2))

# Tests Diagnostics par période
shapiro_pre <- shapiro.test(pre_covid_varx$Res)
bg_pre <- bgtest(pre_covid_varx$Res ~ 1, order = 4) # Breusch-Godfrey au lieu de Ljung-Box (plus robuste)
bp_pre <- bptest(pre_covid_varx$Res ~ pre_covid_varx$Fit)

shapiro_post <- shapiro.test(post_covid_varx$Res)
bg_post <- bgtest(post_covid_varx$Res ~ 1, order = 4)
bp_post <- bptest(post_covid_varx$Res ~ post_covid_varx$Fit)
# Affichage du rapport
cat("\n###############################################################\n")
cat("###        VARX: COMPARISON PRE-COVID vs POST-COVID         ###\n")
cat("###############################################################\n")
cat(sprintf("MAE  Pre-Covid: %.4f  |  MAE  Post-Covid: %.4f\n", mae_pre, mae_post))
cat(sprintf("RMSE Pre-Covid: %.4f  |  RMSE Post-Covid: %.4f\n", rmse_pre, rmse_post))
cat("---------------------------------------------------------------\n")
cat("Normalité (Shapiro) Pre-Covid p-val: ", shapiro_pre$p.value, "\n")
cat("Normalité (Shapiro) Post-Covid p-val:", shapiro_post$p.value, "\n")
cat("Autocorr (Breusch-G) Pre-Covid p-val: ", bg_pre$p.value, "\n")
cat("Autocorr (Breusch-G) Post-Covid p-val:", bg_post$p.value, "\n")
cat("Hétéro (Breusch-P) Pre-Covid p-val:   ", bp_pre$p.value, "\n")
cat("Hétéro (Breusch-P) Post-Covid p-val:  ", bp_post$p.value, "\n")


############################### VECM-X #############################
#!!! pour vecemx: Le VECM est fait pour modéliser des variables Non-Stationnaires qui sont Cointégrées 
#(=elles bougent ensemble sur le long terme).Si vous utilisez le VECM, on repart généralement des variables 
#en Niveau (cad les spreads bruts, pas les diffs).
#!!! et dcp ici je dois reprendre mon df_histo qui a aussi les juste "Spread..." et a ce niveau-ci du code 
#il doit etre en q2 2014-q2 2025 et sans aucun na mais reverifier

# 1. Préparation des données en NIVEAU (Levels)
df_vecmx <- df_histo[, c("Date", 
                         "SpreadGovBY", "SpreadInf", "SpreadUR", 
                         "SpreadGDP", "SpreadTBCommodities", "SpreadCDS", 
                         "SpreadDebt", "SpreadBudgetbalance", "VstoxxEUROPE")]
df_vecmx <- na.omit(df_vecmx)
# 2. Nettoyage des Outliers sur les NIVEAUX (Même logique que pour les diffs) et df_hisot a encore des outliers
# a ce stade
# C'est important car des pics extrêmes peuvent fausser le test de cointégration
df_vecmx <- replace_outliers_iqr(df_vecmx, "SpreadGovBY")
df_vecmx <- replace_outliers_iqr(df_vecmx, "SpreadInf")
df_vecmx <- replace_outliers_iqr(df_vecmx, "SpreadUR")
df_vecmx <- replace_outliers_iqr(df_vecmx, "SpreadGDP")
df_vecmx <- replace_outliers_iqr(df_vecmx, "SpreadTBCommodities")
df_vecmx <- replace_outliers_iqr(df_vecmx, "SpreadCDS")
df_vecmx <- replace_outliers_iqr(df_vecmx, "SpreadDebt")
df_vecmx <- replace_outliers_iqr(df_vecmx, "SpreadBudgetbalance")
df_vecmx <- replace_outliers_iqr(df_vecmx, "VstoxxEUROPE")
# 3. Séparation Endogènes / Exogènes (En Niveau)
# Endogènes : Celles qui interagissent avec tt le systeme: (Spread Gov, CDS)
df_vecmx_endo <- df_vecmx[, c("SpreadGovBY", "SpreadCDS")]
# Exogènes : Les variables macro externes
df_vecmx_exo <- df_vecmx[, c("SpreadGDP", "SpreadUR")]
# 4. Sélection du Lag Optimal
# On utilise VARselect sur les niveaux. Attention, pour un VECM, le lag du VAR sous-jacent est souvent p.
# Le VECM aura alors (p-1) différences.
lag_selection_vecmx <- VARselect(df_vecmx_endo, lag.max = 6, type = "const", exogen = df_vecmx_exo)
opt_lag_vecmx <- lag_selection_vecmx$selection["AIC(n)"]
#Utilise lag.max = 6 : C'est un excellent compromis. Cela permet de capter des effets saisonniers ou cycliques 
#sur un an et demi tout en restant raisonnable pour tes degrés de liberté.
#-> il ma choisi opt_lag_vecmx = 4
cat("\n--- LAG SELECTION (VECMX LEVELS) ---\n")
print(lag_selection_vecmx$selection)
cat(paste("Optimal Lag (AIC):", opt_lag_vecmx, "\n"))
# Sécurité : Johansen a besoin d'au moins K=2 pour fonctionner correctement mathématiquement en VECMX et K etant le 
#opt_lag_vecmx que le code ma proposé dcp ici cest K=4
if (opt_lag_vecmx < 2) opt_lag_vecmx <- 2

# 5. Test de Cointégration de Johansen (Trace Test)
# type = "trace", ecdet = "const" (constante dans la relation de long terme)
# K = opt_lag_vecmx (l'ordre du VAR en niveau)
# dumvar = df_vecmx_exo (les variables exogènes)
johansen_test <- ca.jo(df_vecmx_endo, type = "trace", ecdet = "const", K = opt_lag_vecmx, dumvar = df_vecmx_exo)

summary(johansen_test)
# NOTE IMPORTANTE SUR LE RANG DE COINTÉGRATION (r) :
# Regardez le summary. 
# Si la stat de test > valeur critique (5%) a une ligne, on rejette r=0 /puis r <= 1 /puis ... 
#(donc il y a cointégration a au moins le nombre au dessus a chaque fois).
# Si le test suggère r<=1, vous pouvez changer r=1 ci-dessous.
#je vois que le test me dit quil y a au moins 1 cointegration dcp je fixe rank_r = 1 
rank_r = 1
# 6. Estimation du VECM via OLS (cajorls)
vecm_est <- cajorls(johansen_test, r = rank_r)
summary(vecm_est$rlm) # Pour voir les coefficients détaillés si besoin


# 7. Conversion du VECM en VAR (vec2var) pour le Forecasting
# C'est nécessaire car 'predict' fonctionne mieux sur les objets VAR standards
vecm_to_var <- vec2var(johansen_test, r = rank_r)

# 8. Extraction des Résidus et Fitted Values (Sur le modèle converti)
residuals_vecmx <- as.data.frame(residuals(vecm_to_var))
fitted_vecmx <- as.data.frame(fitted(vecm_to_var))
# Attention aux noms des colonnes dans l'objet converti, ils peuvent changer légèrement et dcp mtn ds residuals_vecmx
# "SpreadGovBY" sappelle mtn "resids of SpreadGovBY" et meme chose pour les 2 autres variables et dcp je dois 
#faire a a chaque fois residuals_vecmx$`resids of SpreadGovBY` pour parler de sa colonne ds residuals_vecmx
# et faire ca fitted_vecmx$`fit of SpreadGovBY`

# --- plot ACTUAL VS FITTED (VECMX) ---
# 1. Création d'un dataframe pour le plot
# On récupère les valeurs ajustées pour la cible
fitted_values_vecmx <- fitted_vecmx$`fit of SpreadGovBY`
# On doit aligner les dates : le vecteur fitted est plus court que le df original
# La longueur est égale au nombre de résidus
n_points <- nrow(fitted_vecmx) 

dates_aligned <- tail(df_vecmx$Date, n_points)
actual_values <- tail(df_vecmx$SpreadGovBY, n_points)
df_plot_vecmx <- data.frame(
  Date = dates_aligned,
  Actual = actual_values,
  Fitted = fitted_values_vecmx
)
# 2. Plot Actual vs Fitted (Style ggplot)
p_vecmx_fit <- ggplot(df_plot_vecmx, aes(x = Date)) +
  # Actual values (Vert foncé continu - pour rappeler le niveau)
  geom_line(aes(y = Actual, color = "Actual Values"), linewidth = 1) +  
  # Fitted values (Rouge pointillé)
  geom_line(aes(y = Fitted, color = "Fitted Values"), linewidth = 1, linetype = "dashed") +  
  labs(
    title = "VECMX: Actual vs. Fitted Values (Spread Gov Bond Yield)",
    subtitle = paste("Model in Levels with Rank r =", rank_r),
    x = "Date",
    y = "Spread (in bps)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  scale_x_yearqtr(format = "%Y Q%q", n = 10) +
  scale_color_manual(values = c("Actual Values" = "darkgreen", "Fitted Values" = "red"))
print(p_vecmx_fit)

# Diagnostics Rapides (Sur la variable cible)
shapiro_vecmx <- shapiro.test(residuals_vecmx$`resids of SpreadGovBY`)
bg_vecmx <- bgtest(residuals_vecmx$`resids of SpreadGovBY` ~ 1, order = 4)
bp_vecmx <- bptest(residuals_vecmx$`resids of SpreadGovBY` ~ fitted_vecmx$`fit of SpreadGovBY`)
mae_vecmx <- mean(abs(residuals_vecmx$`resids of SpreadGovBY`))
rmse_vecmx <- sqrt(mean(residuals_vecmx$`resids of SpreadGovBY`^2))

cat("\n--- DIAGNOSTICS VECMX (Target: SpreadGovBY) ---\n")
cat(sprintf("Shapiro-Wilk (p-val): %.4f\n", shapiro_vecmx$p.value))
cat(sprintf("Breusch-Godfrey (p-val): %.4f\n", bg_vecmx$p.value))
cat(sprintf("Hétéro (Breusch-P) p-val:      %.4f\n", bp_vecmx$p.value))
cat(sprintf("MAE (In Levels): %.4f\n", mae_vecmx))
cat(sprintf("RMSE (In Levels): %.4f\n", rmse_vecmx))



#VECMX FORECASTING (2026-2028)  
# 1. Préparation des Exogènes futures
n_ahead <- 12
# Hypothèse : Last Value Carried Forward (Constante)
last_exo_vecmx <- tail(df_vecmx_exo, 1)
future_exo_vecmx <- last_exo_vecmx[rep(1, n_ahead), ]
rownames(future_exo_vecmx) <- NULL
future_exo_vecmx <- as.matrix(future_exo_vecmx)
# 2. Forecasting avec doubles intervalles (80% et 95%)
# Le modèle vec2var prend 'dumvar' pour les exogènes
fcst_vecmx_95 <- predict(vecm_to_var, n.ahead = n_ahead, dumvar = future_exo_vecmx, ci = 0.95)
#plot(fcst_vecmx_95)
fcst_vecmx_80 <- predict(vecm_to_var, n.ahead = n_ahead, dumvar = future_exo_vecmx, ci = 0.80)

# 3. Création du Dataframe Plot
last_date_hist <- max(df_vecmx$Date)
future_dates <- last_date_hist + seq(0.25, by = 0.25, length.out = n_ahead)

# Extraction pour SpreadGovBY
f_vecmx95 <- fcst_vecmx_95$fcst$SpreadGovBY
f_vecmx80 <- fcst_vecmx_80$fcst$SpreadGovBY
#
df_fcst_vecmx <- data.frame(
  Date = future_dates,
  Forecast = f_vecmx95[, "fcst"],
  Lower95 = f_vecmx95[, "lower"],
  Upper95 = f_vecmx95[, "upper"],
  Lower80 = f_vecmx80[, "lower"],
  Upper80 = f_vecmx80[, "upper"]
)

# 4. Plot VECMX Forecast (En Niveau)
p_vecmx_fcst <- ggplot() +
  # Historique
  geom_line(data = df_vecmx, aes(x = Date, y = SpreadGovBY, color = "Historical (Levels)"), linewidth = 1) +
  # Intervalles
  geom_ribbon(data = df_fcst_vecmx, aes(x = Date, ymin = Lower95, ymax = Upper95, fill = "95% CI"), alpha = 0.3) +
  geom_ribbon(data = df_fcst_vecmx, aes(x = Date, ymin = Lower80, ymax = Upper80, fill = "80% CI"), alpha = 0.5) +
  # Forecast
  geom_line(data = df_fcst_vecmx, aes(x = Date, y = Forecast, color = "Forecast"), linewidth = 1, linetype = "dashed") +
  
  labs(title = "VECMX Forecast: Spread Gov Bond Yield (Levels)",
       subtitle = "Model based on Cointegration (Long-run relationship)",
       y = "Spread (in bps)", x = "Date") +
  scale_color_manual(values = c("Historical (Levels)" = "darkgreen", "Forecast" = "blue")) +
  scale_fill_manual(values = c("95% CI" = "skyblue3", "80% CI" = "steelblue")) +
  theme_gray() +
  theme(legend.position = "bottom") +
  scale_x_yearqtr(format = "%Y Q%q", n = 10)
print(p_vecmx_fcst)

# 5. Comparaison Pre/Post Covid (VECMX)
# Ajustement des dates (vec2var perd opt_lag observations)
# fitted() retourne une matrice, on prend la colonne cible
# Le code ci-dessous ajuste dynamiquement:
n_res <- nrow(residuals_vecmx)
df_eval_vecmx <- data.frame(
  Date = tail(df_vecmx$Date, n_res),
  Res = residuals_vecmx$`resids of SpreadGovBY`
)
pre_cov_vec <- subset(df_eval_vecmx, Date < 2020)
post_cov_vec <- subset(df_eval_vecmx, Date >= 2020)

mae_pre_vec <- mean(abs(pre_cov_vec$Res))
rmse_pre_vec <- sqrt(mean(pre_cov_vec$Res^2))
mae_post_vec <- mean(abs(post_cov_vec$Res))
rmse_post_vec <- sqrt(mean(post_cov_vec$Res^2))

cat("\n###############################################################\n")
cat("###        VECMX: COMPARISON PRE-COVID vs POST-COVID        ###\n")
cat("###############################################################\n")
cat(sprintf("MAE  Pre-Covid: %.4f  |  MAE  Post-Covid: %.4f\n", mae_pre_vec, mae_post_vec))
cat(sprintf("RMSE Pre-Covid: %.4f  |  RMSE Post-Covid: %.4f\n", rmse_pre_vec, rmse_post_vec))
#fin part 3














###############################################################################
###                                                                         ###
###                 Part 4 : Innovation (Kalman Filter / State-Space)       ###
###                                                                         ###
###############################################################################
# 2. PRÉPARATION DES DONNÉES
# Sélection des colonnes (On garde tes variables clés + les Granger Winners)
# On utilise df_histo qui contient déjà toutes tes "diff" calculées plus haut.
#pour kalman nos variables ne doievtn pas dutout etre obligatoirememnt statio mais 
# Création du dataframe final df_kalman pour le Kalman et qui regroupe que les varaibles en diffSpread... et la colonne Date evi
#->mtn df_histo a encore comme en part 1 et part 3 tt les variables et leurs firstdiff en pas cleaned en outliers !!! 
#dcp cest bon je peux enlever les outliers a df_kalman
#et ici on nenleve tjrs pas les outliers de df_histo car cest sur df_kalman quon fait ca

#series in levels
df_kalman <- df_histo[, c("Date", "diffSpreadGovBY", "diffSpreadCDS")]

# Application de la fonction replace_outliers_iqr sur les variables 
df_kalman <- replace_outliers_iqr(df_kalman, "diffSpreadGovBY")
df_kalman <- replace_outliers_iqr(df_kalman, "diffSpreadCDS")
# --- DIAGNOSTIC : VÉRIFICATION DES OUTLIERS REMPLACÉS ---
# 1. Liste des variables à vérifier (toutes sauf la Date)
vars_to_check <- c("diffSpreadGovBY", "diffSpreadCDS")
# 2. Boucle pour comparer df_histo (Avant) et df_kalman (Après)
cat("\n=== RAPPORT DE REMPLACEMENT DES OUTLIERS ===\n")
total_outliers <- 0
for (var in vars_to_check) {
  # On repère les indices où la valeur est différente
  # (On utilise round pour éviter les minuscules erreurs de virgule flottante)
  changed_idx <- which(round(df_histo[[var]], 6) != round(df_kalman[[var]], 6))
  
  if (length(changed_idx) > 0) {
    cat(paste0("\nVariable : ", var, " -> ", length(changed_idx), " modifications.\n"))
    
    # On affiche le détail pour cette variable
    for (i in changed_idx) {
      date_val <- df_kalman$Date[i]
      old_val <- df_histo[[var]][i]
      new_val <- df_kalman[[var]][i]
      
      cat(sprintf("  - %s : Avant = %.4f  ->  Après = %.4f\n", 
                  as.character(date_val), old_val, new_val))
    }
    total_outliers <- total_outliers + length(changed_idx)
  } else {
    cat(paste0("\nVariable : ", var, " -> Aucun outlier détecté.\n"))
  }
}
cat(paste0("\nTotal d'outliers remplacés dans le dataset : ", total_outliers, "\n"))



# =============================================================================
# 3. CONSTRUCTION DU MODÈLE (RÉGRESSION DYNAMIQUE + SAISONNALITÉ)
# =============================================================================
# On prépare la matrice des régresseurs (Exogènes)
X_mat <- as.matrix(df_kalman[, c("diffSpreadCDS")])
# Fonction de construction pour l'optimisation (MLE)
build_kalman_seas <- function(params) {
  # 1. Dynamic Regression + Intercept
  mod_reg <- dlmModReg(X = X_mat, addInt = TRUE)
  # 2. Stochastic Seasonal component
  mod_seas <- dlmModSeas(freq = 4)
  mod <- mod_reg + mod_seas
  
  # 3. Variances (using exp to ensure positivity)
  mod$V <- matrix(exp(params[1]), 1, 1) # Measurement noise
  
  # System Covariance Matrix W
  # params[2]: Intercept drift (Local Level)
  # params[3]: Beta_CDS drift (Time-Varying Parameter)
  # params[4]: Seasonal drift (Allows cycle to evolve)
  mod$W <- diag(c(exp(params[2]), exp(params[3]), exp(params[4]), 0, 0))
  
  # 4. Robust Initialization
  # Using a diffuse prior for C0 ensures the filter isn't biased by starting values
  mod$m0 <- rep(0, 5) 
  mod$C0 <- diag(1e4, 5) # Increased for a more "diffuse" (neutral) start
  return(mod)
}
# =============================================================================
# 4. ESTIMATION DES PARAMÈTRES (MLE)
# =============================================================================
cat("\n--- Estimation des paramètres (Reg Dynamique + Saisonnalité) ---\n")
# 4 paramètres à estimer maintenant (à cause de la saisonnalité)
initial_params <- rep(log(0.01), 4)
fit <- dlmMLE(y = df_kalman$diffSpreadGovBY, parm = initial_params, build = build_kalman_seas)
if(fit$convergence == 0) {
  cat("Convergence successful!\n")
} else {
  cat("Warning: MLE did not fully converge.\n")
}
# Construction du modèle final avec les paramètres optimisés
final_model <- build_kalman_seas(fit$par)
# =============================================================================
# 5. APPLICATION DU FILTRE (FILTERING)
# =============================================================================
filtered <- dlmFilter(df_kalman$diffSpreadGovBY, final_model)
#La fonction dlmFilter calcule le gain de kalman K_t à chaque étape t pour mettre à jour tes 5 états
# Extraction des fitted values
fitted_values <- drop(filtered$f)
# Extracting the latent states (Beta_CDS and Intercept)

# States are in filtered$a (predicted) or filtered$m (filtered)
beta_cds_time_varying <- filtered$m[, 2] # The second column is your Beta_CDS
intercept_time_varying <- filtered$m[, 1] # The first column is your Intercept

# This allows you to plot "How the market sensitive to CDS changed over 10 years"
# =============================================================================
# 6. PRÉPARATION ET FORECASTING (EXTENSION 3 ANS / 12 TRIMESTRES)
# =============================================================================
n_ahead <- 12 
# A. Prévision de la variable exogène (CDS) via ARIMA
fit_cds <- auto.arima(df_kalman$diffSpreadCDS)
fit_cds
plot(df_kalman$diffSpreadCDS)
forecast_cds <- forecast(fit_cds, h = n_ahead)
future_X_values <- as.numeric(forecast_cds$mean)
# Création de la matrice future
future_X <- matrix(future_X_values, nrow = n_ahead, ncol = 1)
colnames(future_X) <- colnames(X_mat)
# B. Combinaison Histoire + Futur pour X
X_full <- rbind(X_mat, future_X)
# C. Combinaison Histoire + Futur pour Y (avec NA pour le futur)
y_full <- c(df_kalman$SpreadGovBY, rep(NA, n_ahead))
# D. Re-construction du modèle pour la prévision
# ATTENTION : Il faut refaire la structure Reg + Seas avec la nouvelle matrice X_full
model_for_forecast <- dlmModReg(X = X_full, addInt = TRUE) + dlmModSeas(freq = 4)
model_for_forecast$V <- final_model$V
model_for_forecast$W <- final_model$W
model_for_forecast$m0 <- rep(0, 5)
model_for_forecast$C0 <- diag(10, 5)

# E. Lancement du Filtre sur la période étendue
filtered_extended <- dlmFilter(y_full, model_for_forecast)
# F. Extraction
all_values <- drop(filtered_extended$f)
fitted_hist <- all_values[1:length(df_kalman$SpreadGovBY)]
forecast_val <- all_values[(length(df_kalman$SpreadGovBY) + 1):length(all_values)]

# 7. VISUALISATION (PLOT RESULTS)
# Dates futures
future_dates <- max(df_kalman$Date) + seq(0.25, by = 0.25, length.out = n_ahead)
# Dataframe Historique
df_plotkalman <- data.frame(
  Date = df_kalman$Date,
  Actual = df_kalman$SpreadGovBY,
  Fitted = as.numeric(fitted_hist)
)
# Dataframe Forecast
df_forecastkalman <- data.frame(
  Date = future_dates,
  Forecast = as.numeric(forecast_val)
)
# Le Graphique
ggplot() +
  # Zone Historique
  geom_line(data = df_plotkalman, aes(x = Date, y = Actual, color = "Actual Values"), linewidth = 1) +
  geom_line(data = df_plotkalman, aes(x = Date, y = Fitted, color = "Fitted Values (Dynamic Reg)"), linewidth = 0.8, linetype = "solid") +
  # Zone Forecast
  geom_line(data = df_forecastkalman, aes(x = Date, y = Forecast, color = "Forecasted Values"), linewidth = 1.2, linetype = "dotted") +
  labs(
    title = "Kalman Filter: Dynamic Regression Model",
    subtitle = "Driven by Time-Varying Sensitivity to CDS",
    x = "Date",
    y = "Change in Spread (bps)",
    color = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c(
    "Actual Values" = "black", 
    "Fitted Values (Dynamic Reg)" = "#E69F00", # Orange
    "Forecasted Values" = "#0072B2" # Bleu
  )) +
  scale_x_yearqtr(format = "%Y Q%q", n = 12)



# 8. DIAGNOSTICS & ÉVALUATION DU MODÈLE KALMAN (PRE VS POST COVID)
# 1. Calcul des Résidus
# On s'assure que le dataframe s'appelle bien df_plotkalman
df_plotkalman$Residuals <- df_plotkalman$Actual - df_plotkalman$Fitted

# 2. Métriques Globales (Toute la période)
mae_global <- mean(abs(df_plotkalman$Residuals), na.rm = TRUE)
rmse_global <- sqrt(mean(df_plotkalman$Residuals^2, na.rm = TRUE))

# Tests Statistiques Globaux
# Normalité
shapiro_test_global <- shapiro.test(df_plotkalman$Residuals)

# Hétéroscédasticité (Breusch-Pagan)
# On teste directement si la variance des résidus dépend des valeurs ajustées (Fitted)
bp_test_global <- bptest(Residuals ~ Fitted, data = df_plotkalman)

# Autocorrélation (Breusch-Godfrey)
# On teste l'autocorrélation pure de la série des résidus (intercept only ~ 1) jusqu'à l'ordre 4
bg_test_global <- bgtest(Residuals ~ 1, order = 4, data = df_plotkalman)

# 3. Séparation Pre-Covid / Post-Covid
pre_covid  <- subset(df_plotkalman, Date < 2020)
post_covid <- subset(df_plotkalman, Date >= 2020)
# 4. Métriques Pre-Covid
mae_pre <- mean(abs(pre_covid$Residuals), na.rm = TRUE)
rmse_pre <- sqrt(mean(pre_covid$Residuals^2, na.rm = TRUE))

shapiro_pre <- shapiro.test(pre_covid$Residuals)
bp_pre <- bptest(Residuals ~ Fitted, data = pre_covid)
bg_pre <- bgtest(Residuals ~ 1, order = 4, data = pre_covid)

# 5. Métriques Post-Covid
mae_post <- mean(abs(post_covid$Residuals), na.rm = TRUE)
rmse_post <- sqrt(mean(post_covid$Residuals^2, na.rm = TRUE))

shapiro_post <- shapiro.test(post_covid$Residuals)
bp_post <- bptest(Residuals ~ Fitted, data = post_covid)
bg_post <- bgtest(Residuals ~ 1, order = 4, data = post_covid)

# =============================================================================
# 9. AFFICHAGE DU RAPPORT DE DIAGNOSTIC
# =============================================================================

cat("\n###############################################################\n")
cat("###             KALMAN FILTER DIAGNOSTIC REPORT             ###\n")
cat("###############################################################\n\n")

cat("--- 1. PERFORMANCE METRICS (ACCURACY) ---\n")
cat(sprintf("Global     : MAE = %.4f | RMSE = %.4f\n", mae_global, rmse_global))
cat(sprintf("Pre-Covid  : MAE = %.4f | RMSE = %.4f\n", mae_pre, rmse_pre))
cat(sprintf("Post-Covid : MAE = %.4f | RMSE = %.4f\n", mae_post, rmse_post))

cat("\n--- 2. NORMALITY OF RESIDUALS (Shapiro-Wilk) ---\n")
cat("H0: Residuals are Normally Distributed (We want p-value > 0.05)\n")
cat(sprintf("Global     : p-value = %.4f %s\n", shapiro_test_global$p.value, ifelse(shapiro_test_global$p.value > 0.05, "(Normal)", "*(Not Normal)*")))
cat(sprintf("Pre-Covid  : p-value = %.4f %s\n", shapiro_pre$p.value, ifelse(shapiro_pre$p.value > 0.05, "(Normal)", "*(Not Normal)*")))
cat(sprintf("Post-Covid : p-value = %.4f %s\n", shapiro_post$p.value, ifelse(shapiro_post$p.value > 0.05, "(Normal)", "*(Not Normal)*")))

cat("\n--- 3. HETEROSCEDASTICITY (Breusch-Pagan) ---\n")
cat("H0: Variance is Constant (Homoscedasticity) (We want p-value > 0.05)\n")
cat(sprintf("Global     : p-value = %.4f %s\n", bp_test_global$p.value, ifelse(bp_test_global$p.value > 0.05, "(Constant Var)", "*(Heteroscedastic)*")))
cat(sprintf("Pre-Covid  : p-value = %.4f %s\n", bp_pre$p.value, ifelse(bp_pre$p.value > 0.05, "(Constant Var)", "*(Heteroscedastic)*")))
cat(sprintf("Post-Covid : p-value = %.4f %s\n", bp_post$p.value, ifelse(bp_post$p.value > 0.05, "(Constant Var)", "*(Heteroscedastic)*")))

cat("\n--- 4. AUTOCORRELATION (Breusch-Godfrey, Order 4) ---\n")
cat("H0: No Serial Correlation up to order 4 (We want p-value > 0.05)\n")
cat(sprintf("Global     : p-value = %.4f %s\n", bg_test_global$p.value, ifelse(bg_test_global$p.value > 0.05, "(No Autocorr)", "*(Autocorrelated)*")))
cat(sprintf("Pre-Covid  : p-value = %.4f %s\n", bg_pre$p.value, ifelse(bg_pre$p.value > 0.05, "(No Autocorr)", "*(Autocorrelated)*")))
cat(sprintf("Post-Covid : p-value = %.4f %s\n", bg_post$p.value, ifelse(bg_post$p.value > 0.05, "(No Autocorr)", "*(Autocorrelated)*")))
cat("\n###############################################################\n")


### 3. INNOVATION ANALYSIS: Time-Varying Parameter (TVP) ###
# Pour prouver l'intérêt du Kalman, on montre que le bêta (sensibilité) change ds le temps
# Extraction des États (Coefficients) lissés
# filtered$m contient l'évolution des coefficients au cours du temps
# Colonne 1 = Intercept, Colonne 2 = Coeff CDS
tvp_coeffs <- drop(filtered$m) 

# Création Dataframe pour le plot (On ignore la première ligne d'initialisation 0)
df_tvp <- data.frame(
  Date = df_kalman$Date,
  Beta_CDS = tvp_coeffs[-1, 2] # Coefficient du CDS
)
# Plot de l'évolution de la sensibilité au CDS
p_beta_cds <- ggplot(df_tvp, aes(x = Date, y = Beta_CDS)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_hline(yintercept = mean(df_tvp$Beta_CDS), linetype="dashed", color="gray") +
  labs(title = "Time-Varying Sensitivity of Spread to CDS",
       subtitle = "Kalman Innovation: Shows how Market Fear (CDS) impact changes over time",
       y = "Beta Coefficient (Impact Value)") +
  theme_minimal() +
  scale_x_yearqtr(format = "%Y Q%q", n = 10)
grid.arrange(p_beta_cds, ncol = 1)




###############################################################################
###           FINAL COMPARISON: ARIMA vs VARX vs KALMAN (INNOVATION)        ###
###############################################################################
# 1. Alignement des longueurs
# Le VARX et le Kalman perdent quelques données au début à cause des Lags/Diffs.
# On doit couper tout le monde à la même longueur pour être juste.
# On prend la longueur la plus courte (celle du VARX généralement)
min_len <- min(length(residualsarima), length(residualsvarx$diffSpreadGovBY), length(df_plotkalman$Residuals))

# On prend les DERNIÈRES observations (les plus récentes) pour aligner cad on part de la derniere et on fait 
#en arriere le min_len de q1 2015-q2 2025
res_arima_clean  <- tail(residualsarima, min_len)
res_varx_clean   <- tail(residualsvarx$diffSpreadGovBY, min_len)
res_kalman_clean <- tail(df_plotkalman$Residuals, min_len)

# 2. Création du Dataframe de Comparaison
# On recalcul MAE/RMSE sur ces vecteurs alignés
metrics_calc <- function(resids) {
  c(MAE = mean(abs(resids)), RMSE = sqrt(mean(resids^2)))
}

perf_arima  <- metrics_calc(res_arima_clean)
perf_varx   <- metrics_calc(res_varx_clean)
perf_kalman <- metrics_calc(res_kalman_clean)

# Tableau Global
df_comparison <- data.frame(
  Model = c("ARIMA (Baseline)", "VARX (Multivariate)", "Kalman Filter (Innovation)"),
  MAE_Global = c(perf_arima["MAE"], perf_varx["MAE"], perf_kalman["MAE"]),
  RMSE_Global = c(perf_arima["RMSE"], perf_varx["RMSE"], perf_kalman["RMSE"])
)
print(df_comparison)

# 3. Comparaison Pre/Post Covid (Sur les vecteurs alignés)
# On doit recréer un vecteur date aligné
#la il prend 42 obs avec mon_len et dcp dates_aligned cest de q1 2015-q2 2025
dates_aligned <- tail(df_part3$Date, min_len)

df_comp_split <- data.frame(
  Date = dates_aligned,
  Res_ARIMA = res_arima_clean,
  Res_VARX = res_varx_clean,
  Res_Kalman = res_kalman_clean
)
# Split
pre_cov_comp <- subset(df_comp_split, Date < 2020)
post_cov_comp <- subset(df_comp_split, Date >= 2020)

# Calcul des RMSE spécifiques
rmse_pre_arima <- sqrt(mean(pre_cov_comp$Res_ARIMA^2))
rmse_pre_varx  <- sqrt(mean(pre_cov_comp$Res_VARX^2))
rmse_pre_kal  <- sqrt(mean(pre_cov_comp$Res_Kalman^2))

rmse_post_arima <- sqrt(mean(post_cov_comp$Res_ARIMA^2))
rmse_post_varx  <- sqrt(mean(post_cov_comp$Res_VARX^2))
rmse_post_kal  <- sqrt(mean(post_cov_comp$Res_Kalman^2))

# Calcul des MAE spécifiques (Mean Absolute Error)
mae_pre_arima <- mean(abs(pre_cov_comp$Res_ARIMA))
mae_pre_varx  <- mean(abs(pre_cov_comp$Res_VARX))
mae_pre_kal   <- mean(abs(pre_cov_comp$Res_Kalman))

mae_post_arima <- mean(abs(post_cov_comp$Res_ARIMA))
mae_post_varx  <- mean(abs(post_cov_comp$Res_VARX))
mae_post_kal   <- mean(abs(post_cov_comp$Res_Kalman))

cat("\n--- BATTLE OF MODELS: RMSE PRE vs POST COVID ---\n")
cat(sprintf("ARIMA  -> Pre: %.4f | Post: %.4f\n", rmse_pre_arima, rmse_post_arima))
cat(sprintf("VARX   -> Pre: %.4f | Post: %.4f\n", rmse_pre_varx, rmse_post_varx))
cat(sprintf("KALMAN -> Pre: %.4f | Post: %.4f\n", rmse_pre_kal, rmse_post_kal))

cat("\n--- BATTLE OF MODELS: MAE PRE vs POST COVID ---\n")
cat(sprintf("ARIMA  -> Pre: %.4f | Post: %.4f\n", mae_pre_arima, mae_post_arima))
cat(sprintf("VARX   -> Pre: %.4f | Post: %.4f\n", mae_pre_varx, mae_post_varx))
cat(sprintf("KALMAN -> Pre: %.4f | Post: %.4f\n", mae_pre_kal, mae_post_kal))


### 2. DIEBOLD-MARIANO TEST (Statistical Significance) ###
cat("\n--- DIEBOLD-MARIANO TESTS (H0: Forecasts have equal accuracy) ---\n")
cat("Alternative H1: Model 2 is less accurate than Model 1 (We want p-value < 0.05)\n\n")

# Test 1: Kalman vs ARIMA
dm_kal_ari <- dm.test(res_kalman_clean, res_arima_clean, alternative = "less", h = 1, power = 2)
cat("1. Kalman vs ARIMA:\n")
print(dm_kal_ari)
if(dm_kal_ari$p.value < 0.05) cat("-> RESULT: Kalman is SIGNIFICANTLY better than ARIMA.\n\n") else cat("-> RESULT: No significant difference.\n\n")

# Test 2: Kalman vs VARX
dm_kal_var <- dm.test(res_kalman_clean, res_varx_clean, alternative = "less", h = 1, power = 2)
cat("2. Kalman vs VARX:\n")
print(dm_kal_var)
if(dm_kal_var$p.value < 0.05) cat("-> RESULT: Kalman is SIGNIFICANTLY better than VARX.\n") else cat("-> RESULT: No significant difference.\n")

# Test 3: VARX vs ARIMA
# On teste si les erreurs du VARX sont "less" (plus petites) que celles de l'ARIMA
dm_var_ari <- dm.test(res_varx_clean, res_arima_clean, alternative = "less", h = 1, power = 2)
cat("3. VARX vs ARIMA:\n")
print(dm_var_ari)
if(dm_var_ari$p.value < 0.05)
  cat("-> RESULT: VARX is SIGNIFICANTLY better than ARIMA.\n") else cat("-> RESULT: No significant difference (VARX implies no improvement over ARIMA).\n")



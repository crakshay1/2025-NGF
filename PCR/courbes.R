# Création de plusieurs dataframes avec nos données de PCR
library(readxl)
library(ggplot2)
library(gridExtra)

# Le chemin est à modifier
chemin <- "C:/Users/taksh/OneDrive - Universite Evry Val d'Essonne/LDD/Biologie/L3/Régulation Expression Génique Eucaryotes/TP/PCR/quantification_PCR.xlsx"

# Un dataframe par gène 
hprt1 <- read_excel(chemin, 
                    col_names = FALSE,
                    sheet = "Feuil1", 
                    range = "A20:K21"
)
hprt1 <- data.frame("Temps" = unlist(hprt1[1,]), "PM" = unlist(hprt1[2,]))
hprt1$Temps <- factor(hprt1$Temps,
                      levels = hprt1$Temps,
                      ordered = TRUE)
hprt1$PM <- as.numeric(hprt1$PM)

atf3 <- read_excel(chemin,
                   col_names = FALSE,
                   sheet = "Feuil1", 
                   range = "A2:F3"
)
atf3 <- data.frame("Temps" = unlist(atf3[1,]), "PM" = unlist(atf3[2,]))
atf3$Temps <- factor(atf3$Temps,
                     levels = atf3$Temps,
                     ordered = TRUE)
atf3$PM <- as.numeric(atf3$PM)

egr1 <- read_excel(chemin, 
                   col_names = FALSE,
                   sheet = "Feuil1", 
                   range = "A8:F9"
)
egr1 <- data.frame("Temps" = unlist(egr1[1,]), "PM" = unlist(egr1[2,]))
egr1$Temps <- factor(egr1$Temps,
                     levels = egr1$Temps,
                     ordered = TRUE)
egr1$PM <- as.numeric(egr1$PM)

cfos <- read_excel(chemin, 
                   col_names = FALSE,
                   sheet = "Feuil1", 
                   range = "A14:F15",
)
cfos <- data.frame("Temps" = unlist(cfos[1,]), "PM" = unlist(cfos[2,]))
cfos$Temps <- factor(cfos$Temps,
                     levels = cfos$Temps,
                     ordered = TRUE)
cfos$PM <- as.numeric(cfos$PM)

# Normalisation via HPRT1
normaliser <- function(gene) {
  newPM <- gene$PM/hprt1$PM[hprt1$Temps %in%  gene$Temps]
  return (newPM)
}

atf3$PM <- normaliser(atf3)
egr1$PM <- normaliser(egr1)
cfos$PM <- normaliser(cfos)

# Résumé 
resume <- data.frame(Temps = hprt1$Temps)
resume <- merge(resume, atf3[, c("Temps", "PM")], by = "Temps", all.x = TRUE)
names(resume)[names(resume) == "PM"] <- "ATF3_PM"
resume <- merge(resume, egr1[, c("Temps", "PM")], by = "Temps", all.x = TRUE)
names(resume)[names(resume) == "PM"] <- "EGR1_PM"
resume <- merge(resume, cfos[, c("Temps", "PM")], by = "Temps", all.x = TRUE)
names(resume)[names(resume) == "PM"] <- "CFOS_PM"

# Tracés des courbes
ggplot(atf3, aes(x = Temps, y = PM)) +
  geom_line(group = 1, color = "purple", size = 1.2) +
  geom_point(color = "purple", size = 3) +
  labs(x = "Temps en minutes", y = "PM") +
  ggtitle("Quantification RT-PCR du gène ATF3 en fonction des temps d'incubation sous NGF") +
  theme_minimal()

ggplot(egr1, aes(x = Temps, y = PM)) +
  geom_line(group = 1, color = "red", size = 1.2) +
  geom_point(color = "red", size = 3) +
  labs(x = "Temps en minutes", y = "PM") +
  ggtitle("Quantification RT-PCR du gène EGR1 en fonction des temps d'incubation sous NGF") +
  theme_minimal()

ggplot(cfos, aes(x = Temps, y = PM)) +
  geom_line(group = 1, color = "darkgreen", size = 1.2) +
  geom_point(color = "darkgreen", size = 3) +
  labs(x = "Temps en minutes", y = "PM") +
  ggtitle("Quantification RT-PCR du gène CFOS en fonction des temps d'incubation sous NGF") +
  theme_minimal()

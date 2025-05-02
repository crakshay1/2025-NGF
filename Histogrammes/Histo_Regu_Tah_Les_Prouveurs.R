makeHistV1 <- function(data){
  # On calcule le Ratio
  data$RatioD <- data$`Nb de Cellules différenciées` / data$`Nb de Cellules totales`
  
  # Moyennes
  plaque3_0_m <- mean(data$RatioD[grepl("p3_0", data$Image)])
  plaque3_20_m <- mean(data$RatioD[grepl("p3_20", data$Image)])
  plaque3_100_m <- mean(data$RatioD[grepl("p3_100", data$Image)])
  plaque4_0_m <- mean(data$RatioD[grepl("p4_0", data$Image)])
  plaque4_20_m <- mean(data$RatioD[grepl("p4_20", data$Image)])
  plaque4_100_m <- mean(data$RatioD[grepl("p4_100", data$Image)])
  
  # Ecart-types
  plaque3_0_e <- sd(data$RatioD[grepl("p3_0", data$Image)])
  plaque3_20_e <- sd(data$RatioD[grepl("p3_20", data$Image)])
  plaque3_100_e <- sd(data$RatioD[grepl("p3_100", data$Image)])
  plaque4_0_e <- sd(data$RatioD[grepl("p4_0", data$Image)])
  plaque4_20_e <- sd(data$RatioD[grepl("p4_20", data$Image)])
  plaque4_100_e <- sd(data$RatioD[grepl("p4_100", data$Image)])
  
  # Création des matrices
  moy_exp <- rbind(
    c(plaque3_0_m, plaque3_20_m, plaque3_100_m),
    c(plaque4_0_m, plaque4_20_m, plaque4_100_m)
  )
  et_exp <- rbind(
    c(plaque3_0_e, plaque3_20_e, plaque3_100_e),
    c(plaque4_0_e, plaque4_20_e, plaque4_100_e)
  )
  
  profilNGF <- c("0 ng/ml", "20 ng/ml", "100 ng/ml")
  rownames(moy_exp) <- c("Plaque 3", "Plaque 4")
  colnames(moy_exp) <- profilNGF
  
  # Plot
  grapheProfil <- barplot(moy_exp, beside = TRUE, names.arg = profilNGF,
                          col = c("#5381E0", "#E05370"),
                          main = "Histogramme des moyennes normalisées\ndes cellules différenciées\npour différentes doses de NGF",
                          ylab = "Moyenne de différenciation",
                          xlab = "Concentration de la dose de NGF",
                          ylim = c(0, 1))
  
  # Ajout des barres d'erreurs
  arrows(grapheProfil, as.vector(moy_exp) - as.vector(et_exp), 
         grapheProfil, as.vector(moy_exp) + as.vector(et_exp), 
         angle = 90, code = 3, length = 0.1)
  
  # Légende
  legend(x = "topright",
         legend = rownames(moy_exp), 
         fill = c("#5381E0", "#E05370"), 
         bty = "n", horiz = TRUE)
}

makeHistV2 <- function(data){
  # On calcule le Ratio
  data$RatioD <- data$`Nb de Cellules différenciées` / data$`Nb de Cellules totales`
  
  # Moyennes
  plaque_0_m <- mean(data$RatioD[grepl("_0", data$Image)])
  plaque_20_m <- mean(data$RatioD[grepl("_20", data$Image)])
  plaque_100_m <- mean(data$RatioD[grepl("_100", data$Image)])
  
  # Ecart-types
  plaque_0_e <- sd(data$RatioD[grepl("_0", data$Image)])
  plaque_20_e <- sd(data$RatioD[grepl("_20", data$Image)])
  plaque_100_e <- sd(data$RatioD[grepl("_100", data$Image)])
  
  # Création des matrices
  moy_exp <- rbind(
    c(plaque_0_m, plaque_20_m, plaque_100_m)
  )
  et_exp <- rbind(
    c(plaque_0_e, plaque_20_e, plaque_100_e)
  )
  
  profilNGF <- c("0 ng/ml", "20 ng/ml", "100 ng/ml")
  rownames(moy_exp) <- c("Moyenne des plaques")
  colnames(moy_exp) <- profilNGF
  
  # Plot
  grapheProfil <- barplot(moy_exp, beside = TRUE, names.arg = profilNGF,
                          col = "#E05370",
                          main = "Histogramme des moyennes normalisées\ndes cellules différenciées\npour différentes doses de NGF",
                          ylab = "Moyenne de différenciation",
                          xlab = "Concentration de la dose de NGF",
                          ylim = c(0, 1))
  
  # Ajout des barres d'erreurs
  arrows(grapheProfil, as.vector(moy_exp) - as.vector(et_exp), 
         grapheProfil, as.vector(moy_exp) + as.vector(et_exp), 
         angle = 90, code = 3, length = 0.1)
  
  # Légende
  legend(x = "topright",
         legend = rownames(moy_exp), 
         fill = "#E05370", 
         bty = "n", horiz = TRUE)
}

makeWholeHist <- function(data) {
  # Pour voir l'ensemble de l'expérience
  data$RatioD <- data$`Nb de Cellules différenciées` / data$`Nb de Cellules totales`
  data$RatioND <- 1-data$RatioD
  moyennes <- c()
  ecart_types <- c()
  moyennes1 <- c()
  ecart_types1 <- c()
  for (i in seq(from = 1, to = length(data$`Nom  Binôme`)+5, by = 5)) {
    moyennes <- c(moyennes, mean(data[i:(i+4),]$RatioD))
    ecart_types <- c(ecart_types, sd(data[i:(i+4),]$RatioD))
    moyennes1 <- c(moyennes1, mean(data[i:(i+4),]$RatioND))
    ecart_types1 <- c(ecart_types1, sd(data[i:(i+4),]$RatioND))
  }
  moyennes <- moyennes[-13]
  ecart_types <- ecart_types[-13]
  moyennes1 <- moyennes1[-13]
  ecart_types1 <- ecart_types1[-13]
  
  MOY <- rbind(moyennes, moyennes1)
  ET <- rbind(ecart_types,ecart_types1)
  rownames(MOY) <- c("Différenciées", "Non Différenciées")
  colnames(MOY) <- c("p3_0", "p3_20", "p3_100", "p3_0_1",
                     "p3_20_1", "p3_100_1", "p4_0", "p4_20",
                     "p4_100", "p4_0_1", "p4_20_1", "p4_100_1"
  )
  grapheProfils <- barplot(MOY, beside =TRUE, names.arg = colnames(MOY),
                           col = c("#d67b16", "#01538f"), #Les couleurs de la fac RAAAH
                           main = "Histogramme des moyennes normalisées des cellules", 
                           ylab = "Moyenne",
                           xlab = "Conditions auxquelles les cellules sont soumises",
                           ylim = c(0,(max(MOY))+0.1))
  arrows(grapheProfils, MOY - ET, grapheProfils,MOY + ET, 
         angle= 90, code =3, length =0.01)
  par(mar = c(6, 4, 4, 2) + 0.1)
  legend("bottomleft", inset = c(0.05, -0.225), xpd = TRUE,
         legend = rownames(MOY), 
         fill = c("#d67b16", "#01538f"), 
         bty = "n", 
         horiz = TRUE)
}

###########################################################################

# Création d'un dataframe avec toutes les données (MERCI TASSILIE)
library(readxl)
# Le chemin est à modifier
chemin <- "C:/Users/taksh/OneDrive - Universite Evry Val d'Essonne/LDD/Biologie/L3/Régulation Expression Génique Eucaryotes/TP/Photos/Photos.xlsx"
exps <- read_excel(chemin, 
                   sheet = "Comptage", 
                   range = "A1:D61"
)
# Génération de l'histogramme
makeHistV1(exps) # Figure inutilisée
makeHistV2(exps) # Figure 2
makeWholeHist(exps) # Figure inutilisée

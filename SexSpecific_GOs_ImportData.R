# Cargar las librerías dplyr, tidyr y ggplot2
library(dplyr)
library(tidyr)
library(ggplot2)

# Establecer el directorio de trabajo
setwd("./../Data/GO_enrichment/")

# IMPORTACIÓN DE DATOS
###############################################################
# Se evaluaron todos los TFs presentes en los promotores de DMC1, DUET, SPO11-1, UBC22 y FST.
# Se almacenaron todos los términos GOs de PCBase 2.0 en archivos .txt para cada TF.
data <- read.csv("./JAG_GO.txt", header = FALSE, sep = "\t") # Abrir archivo .txt
head(data) # Revisar las etiquetas de columnas

# Dentro de PCBase 2.0, los términos GO vienen acompañados de diferentes columnas con las siguientes combinaciones.
# Se emplearon las combinaciones adecuadas para nombrar cada columna.
colnames(data) <- c("GO_ID", "GO_category", "GO_term", "All", "JA", "MeJA", "Nontreatment", "Others")
colnames(data) <- c("GO_ID", "GO_category", "GO_term", "All", "FR irradiation")
colnames(data) <- c("GO_ID", "GO_category", "GO_term", "All", "DSG", "Nontreatment")
colnames(data) <- c("GO_ID", "GO_category", "GO_term", "All", "FR_irradiation")
colnames(data) <- c("GO_ID", "GO_category", "GO_term", "All", "ABA", "Dark")
colnames(data) <- c("GO_ID", "GO_category", "GO_term", "All", "Nontreatment")

# Filtrar GO IDs de meiosis femenina y masculina
GO_meiotic_numbers <- c("GO:0007147", "GO:0007142") # Crear vector con GO IDs para ambos sexos
data_filtered <- data %>% filter(GO_ID %in% GO_meiotic_numbers)
# All incluye los valores de significancia de GOs para todas las condiciones experimentales
data_filtered <- data_filtered[c("GO_ID", "GO_category", "GO_term", "All")]

TF <- "JAG" # Nombrar el TF bajo análisis
# Añadir una columna con el nombre del TF
data_filtered$TF_name <- c(rep(TF, nrow(data_filtered)))
data_filtered

# Recuperar data_filtered en un archivo .txt
# El nombre se autocompleta según el valor de TF
path <- file.path("./../../Results/TF_meioticGOs/", paste(TF, "_meioticGOs.txt", sep = ""))
write.csv(data_filtered, file = path, row.names = FALSE)
# Cargar las librerías dplyr, tidyr y ggplot2
library(dplyr)
library(tidyr)
library(ggplot2)

# fijar directorio de trabajo
setwd(dir = "C:/Users/user/Desktop/MASTER_PROJECT/AVANCES/Bioinformatic_analysis/Data_visualization/Promoter Analysis - R/PCDBase/Data/")

# Importar datos obtenidos de PCBase 2.0
pDMC1 <- read.csv(file = "./promoters/pDMC1_500toTSS.txt", header = TRUE, sep = "\t")
pSPO11 <- read.csv(file = "./promoters/pSPO11-1_500toTSS.txt", header = TRUE, sep = "\t")
pDUET <- read.csv(file = "./promoters/pDUET_500toTSS.txt", header = TRUE, sep = "\t")
pUBC22 <- read.csv(file = "./promoters/pUBC22_500toTSS.txt", header = TRUE, sep = "\t")
pFST <- read.csv(file = "./promoters/pFST_500toTSS.txt", header = TRUE, sep = "\t")

# Añadir una columna con los nombres de cada promotor
pDMC1$promoter <- c(rep("pDMC1", nrow(pDMC1)))
pSPO11$promoter <- c(rep("pSPO11-1", nrow(pSPO11)))
pDUET$promoter <- c(rep("pDUET", nrow(pDUET)))
pUBC22$promoter <- c(rep("pUBC22", nrow(pUBC22)))
pFST$promoter <- c(rep("pFST", nrow(pFST)))
head(pDMC1) # verificar el dataframe pDMC1 con todas las columnas

# Estimar las frecuencias de los TFs que se unen a cada promotor
freq_data <- pDMC1 %>%
  count(`Protein.ID`, `promoter`) %>%
  rename(Frequency = n) %>%
  arrange(desc(Frequency))  # Reordenar los TFs de mayor a menor frecuencia
freq_data # Verificar el dataframe

# Recuperar los datos en un archivo .txt
write.csv(x = freq_data, file = "./../Results/15_mostcommonTFsbyPromoter/pDMC1_AllTFsID.txt", row.names = FALSE)


# Filtrar los 15 TFs más abundantes según promotor
# Crear la funcion top15
top15 <- function(datos){
  freq_data <- datos %>%
    count(`Protein.ID`, `promoter`) %>%
    rename(Frequency = n) %>%
    arrange(desc(Frequency)) %>%
    head(15)
  return(freq_data)
}

# Crear un data frame con los 15 TFs más abundantes en todos los promotores
alldata <- rbind(top15(pDMC1), top15(pSPO11), top15(pDUET), top15(pUBC22), top15(pFST))
alldata
# Recuperar alldata en un archivo .txt
write.csv(x = alldata, file = "./../Results/15_mostcommonTFsbyPromoter/15mostCommonIn5promoters.txt", row.names = FALSE)

# Agrupar por nombre de TF, sumar sus respectivas frecuencias y seleccionar los 15 TFs más abundantes
freq_alldata <- alldata %>%
  group_by(Protein.ID) %>%
  summarise(frecuencia_total = sum(Frequency)) %>%
  arrange(desc(frecuencia_total)) %>%  # Ordenar de mayor a menor
  top_n(15, frecuencia_total) %>%
  ungroup()  # Eliminar el agrupamiento para no interferir en otras operaciones
freq_alldata

# Crear un vector con los 15 TFs más abundantes en freq_alldata
top_15TFs <- as.vector(freq_alldata$Protein.ID)

# Seleccionar los TFs de alldata que se encuentran dentro de top15TFs
data_for_plot <- alldata[alldata$Protein.ID %in% top_15TFs,]

# Añadir orden a data_for_plot según el orden mostrado en top15TFs
data_for_plot$Protein.ID <- factor(data_for_plot$Protein.ID, levels = top_15TFs)

# Crear un vector de colores personalizados para cada promotor
custom_colors <- c("darkred", "#00A86B", "#62CDFF", "#BD1A1A", "#EDA6C4")

# Crear el gráfico de barras apiladas
ggplot(data_for_plot, aes(x = Protein.ID, y = Frequency, fill = promoter)) +
  geom_bar(stat = "identity", position = "stack") +  # Apilar las barras
  scale_x_discrete( labels = c("HY5", "MYB3R1", "TREE1", "MYB3R4", "WRKYs", "PIF4", "MYC2", "ARR1", "NFYB2", "HB7", "NFYC2", "bZIP68", "AT5G04760", "ANAC032", "REF6")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Transcription factor", y = "TFs number", fill = "Promoter") +  # Etiquetas
  theme_minimal() +  # Tema minimalista
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, face = "bold", colour = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank()
    )








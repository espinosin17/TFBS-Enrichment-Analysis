# Fijar directorio de trabajo
setwd(dir = "C:/Users/user/Desktop/MASTER_PROJECT/AVANCES/Bioinformatic_analysis/Data_visualization/Promoter Analysis - R/PCDBase/Results/")

# Importar datos de OG para los 29 TFs meióticos en un dataframe
# Abrir todos los archivos con los GOs de ambos sexos ubicados en el directorio TF_meioticGOs
meioticTFsList <-  lapply(list.files("./TF_meioticGOs/", pattern="*.txt", full.names=TRUE), read.csv)%>%
  lapply(function(df){
    colnames(df)[4] <- "-log(p-value)" # Cambiar el nombre de la columna All
    return(df)
  }) %>%
  do.call(rbind, .) # Unir todas las observaciones en un solo dataframe
unique(meioticTFsList$TF_name) # Corroborar los 29 TFs

data <- meioticTFsList #Crear otro objeto con los mismas observaciones

################################################################################
## ANÁLISIS DE LA MEIOSIS FEMENINA
# Filtrar los TFs con funciones putativas en la meiosis femenina
data_filtered <- filter(data, GO_ID == "GO:0009554")

# Crear gráfico de barras horizontales
# Ordenar el eje x en función de -log(p-value)
ggplot(data_filtered, aes(x = reorder(`TF_name`, `-log(p-value)`), y = `-log(p-value)`)) +
  geom_col(fill = "#90EE90") +
  scale_x_discrete( labels = c("DREB2A", "ABF1", "REF6", "ARR1", "AZF1", "NFYC2", "GBF2",
                               "PIF1", "SPT6L", "LEC1", "ABF3", "ANAC032", "RD26", "MYB44", 
                               "SEP3", "JAG")) +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", size = 0.8) +  # Línea horizontal
  coord_flip() +  # Invertir los ejes x e y
  scale_y_continuous(breaks = seq(0, max(data_filtered$`-log(p-value)`), by = 0.5), expand = c(0,0)) +
  labs(title = "Female meiosis (GO:0009554)", 
       x = "Transcription factors", 
       y = "-log(p-value)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, face = "bold", colour = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  )

################################################################################
## ANÁLISIS DE LA MEIOSIS MASCULINA
# Filtrar los TFs con funciones putativas en la meiosis masculina
data_filtered <- filter(data, GO_ID == "GO:0009556")

# Definir un nivel de significancia alpha = 0.05
# Los TFs con significancia superior a 0.05 recibirán la etiqueta "superior" y viceversa
data_filtered$color_group <- ifelse(data_filtered$`-log(p-value)` > -log10(0.05), "superior", "inferior")

# Graficar
ggplot(data_filtered, aes(x = reorder(`TF_name`, `-log(p-value)`), y = `-log(p-value)`)) +
  geom_col(fill = "#90EE90") +
  scale_x_discrete( labels = c("GBF2", "SEP3", "LEC1", "FHY1", "MYB3R1", "PhyA")) +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed", size = 0.8) +  # Línea horizontal
  coord_flip() +  # Invierte los ejes
  scale_y_continuous(breaks = seq(0, max(data_filtered$`-log(p-value)`), by = 0.5), expand = c(0,0)) +
  labs(title = "Male meiosis (GO:0009556)", 
       x = "Transcription factors", 
       y = "-log(p-value)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", colour = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, face = "bold", colour = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  )
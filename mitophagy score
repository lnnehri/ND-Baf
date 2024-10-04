# Yeni veri setine göre df oluşturma
data <- data.frame(
  group = c(2, 3, 4),
  
  PINK1 = c(1.168312444, 1.209997078, 1.641245948), # PINK1 inducer doi: 10.1242/jcs.093849
  SRC = c(1.170263989, 0.863352332, 1.238465989), # SRC supressor https://doi.org/10.3389/fphar.2022.897046
  ULK1 = c(1.89798269, 1.552060398, 2.596375092), # ULK1 inducer https://doi.org/10.1038/s41598-021-00170-4
  MTERF3 = c(-1.037074875, NA, -1.2189022), # MTERF3 inducer https://doi.org/10.3389/fonc.2023.1132559
  PGAM5 = c(-1.062996008, NA, NA), # PGAM5 #inducer https://doi.org/10.1016/j.toxlet.2022.10.003
  TOMM40 = c(-1.314918505, NA, -1.268570071), # TOMM40 #inducer doi: https://doi.org/10.1101/2023.06.24.546411
  TOMM5 = c(-1.564804791, NA, -1.2770355), # TOMM5 inducer https://www.genecards.org/cgi-bin/carddisp.pl?gene=TOMM5
  TOMM6 = c(-1.205514409, NA, -0.8967824), # TOMM6 inducer https://www.uniprot.org/uniprotkb/Q9CQN3/entry
  UBE2N = c(-1.608133738, NA, -1.6661424), # UBE2N inducer  DOI: 10.1242/jcs.146035
  VDAC1 = c(-1.002400677, NA, -1.4883174), # VDAC1 inducer https://doi.org/10.1038/ncb2012
  MAP1LC3B = c(NA, NA, 1.16986854), # MAP1LC3B inducer DOI: 10.1080/15548627.2022.2077552
  OPTN = c(0.763209199, 1.1236032, 1.146394842), # OPTN inducer doi: 10.1080/15548627.2015.1009792.
  SQSTM1 = c(NA, 0.94590932, 1.467901722), # SQSTM1 inducer  doi: 10.1080/15548627.2019.1643185.
  FUNDC1 = c(NA, NA, -1.1838776) # FUNDC1 #inducer https://doi.org/10.3389/fphar.2022.897046
) 

# NA değerleri 0 ile değiştirme
data[is.na(data)] <- 0


# normalization
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

normalize1 <- function(x) {
  return(1 - ((x - min(x)) / (max(x) - min(x))))
}


#supressors
data$normalized_SRC <- normalize1(data$SRC) 
#inducers
data$normalized_PINK1 <- normalize(data$PINK1) 
data$normalized_ULK1 <- normalize(data$ULK1) 
data$normalized_MTERF3 <- normalize(data$MTERF3) 
data$normalized_PGAM5 <- normalize(data$PGAM5) 
data$normalized_TOMM40 <- normalize(data$TOMM40) 
data$normalized_TOMM5 <- normalize(data$TOMM5) 
data$normalized_TOMM6 <- normalize(data$TOMM6) 
data$normalized_UBE2N <- normalize(data$UBE2N) 
data$normalized_VDAC1 <- normalize(data$VDAC1) 
data$normalized_MAP1LC3B <- normalize(data$MAP1LC3B) 
data$normalized_OPTN <- normalize(data$OPTN) 
data$normalized_SQSTM1 <- normalize(data$SQSTM1) 
data$normalized_FUNDC1 <- normalize(data$FUNDC1) 

data$Adjusted_Score <- data$normalized_PINK1 + data$normalized_SRC + 
  data$normalized_ULK1 + data$normalized_MTERF3 + data$normalized_PGAM5 + 
  data$normalized_TOMM40 + data$normalized_TOMM5 + data$normalized_TOMM6 + 
  data$normalized_UBE2N + data$normalized_VDAC1 + data$normalized_MAP1LC3B + 
  data$normalized_OPTN + data$normalized_SQSTM1 + data$normalized_FUNDC1



# ROS ayrimi
data$category <- ifelse(data$Adjusted_Score > mean(data$Adjusted_Score), "more mitophagy", "less mitophagy")
print(data)

View(data)

for (i in 1:nrow(data)) {
  cat("Group:", data$group[i], "\t", 
      "Category:", data$category[i], "\n")
}

library('ggplot2')
ggplot(data, aes(x = group, y = Adjusted_Score, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("more mitophagy" = "red", "less mitophagy" = "blue")) +
  labs(x = "Group", y = "Adjusted Score", fill = "Category") +
  theme_minimal()

#necrosis score
# df olusturma
data <- data.frame(
  group = c(2, 3, 4),
  ARHGEF2 = c(1.321720956, NA, 1.974636925), # inducer https://www.uniprot.org/uniprotkb/Q60875/entry
  BIRC3 = c(-1.291467683, NA, -1.077535467), #inducer https://doi.org/10.1038/srep21710
  MLKL = c(-1.141945718, NA, -1.01753728), #inducer doi: 10.1093/jmcb/mjaa055
  RBCK1 = c(NA, 1.004551771, 1.147692203), #supressor DOI: 10.1074/jbc.M701913200
  YBX3 = c(-1.24827491, NA, -1.061638287) #supressor https://www.ncbi.nlm.nih.gov/gene/8531
  
)

data[is.na(data)] <- 0
# bakiyorum
print(data)


# normalization
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

normalize1 <- function(x) {
  return(1 - ((x - min(x)) / (max(x) - min(x))))
}


# inducers ornek:
# data$normalized_NOS3 <- normalize(data$NOS3) # inducer doi: 10.1172/JCI21968

#supressors ornek:
# data$normalized_WASL <- normalize1(data$WASL) # suppressor doi: 10.1016/j.parkreldis.2021.02.001

#asagidaki kodu incuder mi suppressor mu diye duzenle, supressor icin normalize1 yap fonksiyonu
data$normalized_ARHGEF2 <- normalize(data$ARHGEF2)

data$normalized_BIRC3 <- normalize(data$BIRC3)

data$normalized_MLKL <- normalize(data$MLKL)

data$normalized_RBCK1 <- normalize1(data$RBCK1)

data$normalized_YBX3<- normalize1(data$YBX3)


###
#adjusted score creation
data$Adjusted_Score <- data$normalized_ARHGEF2 + data$normalized_BIRC3 + 
  data$normalized_MLKL + data$normalized_RBCK1 + data$normalized_YBX3

# ROS ayrimi
data$category <- ifelse(data$Adjusted_Score > mean(data$Adjusted_Score), "more necrosis", "less necrosis")
print(data)

View(data)

for (i in 1:nrow(data)) {
  cat("Group:", data$group[i], "\t", 
      "Category:", data$category[i], "\n")
}

library('ggplot2')
ggplot(data, aes(x = group, y = Adjusted_Score, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("more necrosis" = "red", "less necrosis" = "blue")) +
  labs(x = "Group", y = "Adjusted Score", fill = "Category") +
  theme_minimal()

#pyroptosis score
# df olusturma
#, , , , , , , 


data <- data.frame(
  group = c(2, 3, 4),
  
  TP53 = c(NA, NA, 1.181459494), #inducer  doi: 10.1155/2019/8746895
  
  CASP5 = c(2.430418128, NA, 2.962821488), #inducer https://doi.org/10.1186/s12890-023-02456-x
  
  CYCS = c(-1.542602157, NA, NA), #supressor doi: 10.3389/fcell.2020.630771
  
  IRF2 = c(1.362099859, NA, 1.102453897), #inducer doi: 10.1126/scisignal.aax4917.
  
  CASP1 = c(3.381815532, NA, NA), #inducer https://doi.org/10.1186/s12890-023-02456-x
  
  IL1A = c(NA, NA, 1.841566059), #inducer  doi: 10.1016/j.isci.2020.101070
  
  IL18 = c(1.114598545, 0.970548107, NA), #inducer  doi: 10.1371/journal.ppat.1002452
  
  HMGB1 = c(-1.498326123, -1.053874525, -2.299779688) #inducer https://doi.org/10.1038/s41467-020-18443-3

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
data$normalized_TP53<- normalize(data$TP53)

data$normalized_CASP5<- normalize(data$CASP5)

data$normalized_CYCS <- normalize1(data$CYCS)

data$normalized_IRF2 <- normalize(data$IRF2)

data$normalized_CASP1 <- normalize(data$CASP1)

data$normalized_IL1A <- normalize(data$IL1A)

data$normalized_HMGB1 <- normalize(data$HMGB1)

data$normalized_IL18 <- normalize(data$IL18)


###
#adjusted score creation
data$Adjusted_Score <- data$normalized_TP53 + data$normalized_CASP5 + 
  data$normalized_CYCS + data$normalized_IRF2 + data$normalized_CASP1 + 
  data$normalized_IL1A+ data$normalized_HMGB1 + data$normalized_IL18 
# pyroptosis ayrimi
data$category <- ifelse(data$Adjusted_Score > mean(data$Adjusted_Score), "more pyroptosis", "less pyroptosis")
print(data)

View(data)

for (i in 1:nrow(data)) {
  cat("Group:", data$group[i], "\t", 
      "Category:", data$category[i], "\n")
}

library('ggplot2')
ggplot(data, aes(x = group, y = Adjusted_Score, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("more pyroptosis" = "red", "less pyroptosis" = "blue")) +
  labs(x = "Group", y = "Adjusted Score", fill = "Category") +
  theme_minimal()

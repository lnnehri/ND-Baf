#ferroptosis score
# df olusturma
data <- data.frame(
  group = c(2, 3, 4),
  AKR1C1 = c(2.936718154, 2.673912804, 3.738079368), #suppressor https://doi.org/10.1038/s41422-020-00441-1
  
  AKR1C2 = c(2.995378094, 2.046183473, 3.694122759), #supressor https://doi.org/10.1038/s41422-020-00441-1
  
  AKR1C3 = c(1.985725653, 1.244414976, 1.842055813), #supressor https://doi.org/10.1038/s41422-020-00441-1
  
  CTH = c(1.20123897, 0.890424175, 2.052880263), #inducer https://doi.org/10.1016/j.ymthe.2021.03.022
  
  GCLC = c(NA, NA, 1.356673486), #supressor doi: 10.1016/j.cmet.2020.12.007.
  
  HMOX1 = c(NA, NA, 1.206506022), #inducer DOI: 10.1080/02713683.2022.2138450
  
  MAP1LC3B = c(NA, NA, 1.169868537), #inducer https://doi.org/10.1038/s41598-024-54837-9
  
  SAT1 = c(0.879203951, 0.8981564, 1.221946723), # inducer DOI: 10.1073/pnas.1607152113
  
  SAT2 = c(1.586680733, NA, 1.95748785), # inducer https://doi.org/10.1038/s41420-023-01371-8
  
  SLC3A2 = c(1.435279674, NA, 1.725888771), #supressor doi: 10.1186/s41065-022-00225-0
  
  SLC7A11 = c(0.950895624, NA, 1.26086611), #supressor https://doi.org/10.1038/s41467-023-39401-9
  
  TP53 = c(NA, NA, 1.181459494), # inducer https://doi.org/10.1038/s41420-023-01517-8
  
  GSS = c(-0.803089308, NA, -1.217891148), # suppressor https://doi.org/10.1038/s41419-023-06359-x
  
  HSPB1 = c(-1.354869665, NA, -1.358459408), #supressor https://doi.org/10.1038/s41419-023-05972-0
  
  LPCAT3 = c(NA, NA, -1.071017512), #inducer https://doi.org/10.1038/s41392-020-00216-5
  
  NOX1 = c(-2.98504751, -1.921708293, -4.35222649), #inducer https://doi.org/10.1038/s41422-020-00441-1
  
  ALOX15 = c(1.798662566, NA, NA), #inducer https://doi.org/10.1038/s41392-022-01090-z
  
  ACSL4 = c(-1.29656816, NA, NA), #inducer https://doi.org/10.1038/nchembio.2239
  
  FDFT1 = c(-1.105247112, NA, NA), #inducer  DOI: 10.1002/cam4.4716
  
  GCLM = c(-1.178905203, NA, NA), #supressor  doi: 10.1074/jbc.RA119.009548
  
  ACSL5 = c(NA, 1.132628004, NA) #supressor https://doi.org/10.1038/s41392-024-01769-5
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
data$normalized_CTH <- normalize(data$CTH) 
data$normalized_HMOX1 <- normalize(data$HMOX1) 
data$normalized_MAP1LC3B <- normalize(data$MAP1LC3B) 
data$normalized_SAT1 <- normalize(data$SAT1) 
data$normalized_SAT2 <- normalize(data$SAT2) 
data$normalized_TP53 <- normalize(data$TP53) 
data$normalized_LPCAT3<- normalize(data$LPCAT3) 
data$normalized_NOX1<- normalize(data$NOX1) 
data$normalized_ALOX15<- normalize(data$ALOX15) 
data$normalized_ACSL4<- normalize(data$ACSL4) 
data$normalized_FDFT1<- normalize(data$FDFT1) 

#supressors ornek:
data$normalized_AKR1C1 <- normalize1(data$AKR1C1) 
data$normalized_AKR1C2 <- normalize1(data$AKR1C2) 
data$normalized_AKR1C3 <- normalize1(data$AKR1C3) 
data$normalized_GCLC <- normalize1(data$GCLC)  
data$normalized_SLC3A2 <- normalize1(data$SLC3A2)  
data$normalized_SLC7A11<- normalize1(data$SLC7A11) 
data$normalized_GSS<- normalize1(data$GSS) 
data$normalized_HSPB1<- normalize1(data$HSPB1) 
data$normalized_GCLM<- normalize1(data$GCLM) 
data$normalized_ACSL5<- normalize1(data$ACSL5) 



###
#adjusted score creation
data$Adjusted_Score <- data$normalized_AKR1C1 + data$normalized_AKR1C2 + 
  data$normalized_AKR1C3 + data$normalized_CTH + data$normalized_GCLC + 
  data$normalized_HMOX1 + data$normalized_MAP1LC3B + data$normalized_SAT1 + 
  data$normalized_SAT2 + data$normalized_SLC3A2 + data$normalized_SLC7A11 + 
  data$normalized_TP53 + data$normalized_GSS + data$normalized_HSPB1 + 
  data$normalized_LPCAT3 + data$normalized_NOX1 + data$normalized_ALOX15 + 
  data$normalized_ACSL4 + data$normalized_FDFT1 + data$normalized_GCLM + 
  data$normalized_ACSL5

# ROS ayrimi
data$category <- ifelse(data$Adjusted_Score > mean(data$Adjusted_Score), "more ferrotopsis", "less ferrotopsis")
print(data)

View(data)

for (i in 1:nrow(data)) {
  cat("Group:", data$group[i], "\t", 
      "Category:", data$category[i], "\n")
}

library('ggplot2')
ggplot(data, aes(x = group, y = Adjusted_Score, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("more ferrotopsis" = "red", "less ferrotopsis" = "blue")) +
  labs(x = "Group", y = "Adjusted Score", fill = "Category") +
  theme_minimal()

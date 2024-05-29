#R code for apoptosis scoring
#NR as control
#as groups: 2 for ND, 3 for BAF, 4 for ND+BAF
#inhibitors as - and inducers are +
#ENDOD1 gene is eliminated due to un_deterministic effect in terms of inhibition or inducing
#logFC2 values are used to define gene lists
#only threshold is p_adj<0.01 genes, no threshold for logFC2 values

data <- data.frame(
  group = c(2, 3, 4),
  #inducer BCL2L11, 10018
  BCL2L11= c(1.64245599172252, 0.944772245547861, 1.71307243782403),
  #inducer TP53INP1, 94241
  TP53INP1 = c(1.61370266713904, 1.41775168183201, 1.74967782125782), 
  #supressor BIRC3, 330
  BIRC3 = c(-1.29146768268436, 0, -1.07753546657281),
  #inducer CYCS, 54205
  CYCS = c(-1.54260215746871, 0, 0),
# ENDOD1 = c(), 23052
  #supressor TNFRSF10D, 8793
  TNFRSF10D = c(-1.04784716269947, 0, 0),
  #supressor PIK3CD, 5293
  PIK3CD = c(0, 1.41383716088679, 0),
  #supressor PIK3R3, 8503
  PIK3R3 = c(0, 0, 1.91423388178196),
  #supressor BCL2, 596
  BCL2 = c(0, 0, -1.49750385118025),
  #supressor BCL2L1, 598
  BCL2L1 = c(0, 0, -1.61987175055399),
  #supressor CAPN2, 824
  CAPN2 = c(0, 0, -1.19188452024524),
  #inducer PPP3CB, 5532
  PPP3CB = c(0, 0, -1.59985351463316),
  #inducer PRKX, 5613
  PRKX = c(-0.716682839874494, 0, -1.73740842525519),
  #inducer ATM, 472
  ATM  = c(0, 0, 1.01540458680992),
  #inducer IL1A, 3552
  IL1A  = c(0, 0, 1.84156605899359),
  #inducer IRAK2, 3656
  IRAK2  = c(0, 0, 1.11311074859038),
  #inducer PRKACB, 5567
  PRKACB  = c(0, 0, 1.28136556618915),
  #inducer TP53, 7157
  TP53  = c(0, 0, 1.1814594944741),
  #supressor AKT1, 207
  AKT1 = c(0, 0, -1.12537252983054)
)

# normalization for inducers
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# normalization for supressors
normalize1 <- function(x) {
  return(1 - ((x - min(x)) / (max(x) - min(x))))
}

#for inducers

data$normalized_BCL2L11 <- normalize(data$BCL2L11)
data$normalized_TP53INP1 <- normalize(data$TP53INP1)
data$normalized_CYCS <- normalize(data$CYCS)
data$normalized_PPP3CB <- normalize(data$PPP3CB)
data$normalized_PRKX <- normalize(data$PRKX)
data$normalized_ATM <- normalize(data$ATM)
data$normalized_IL1A <- normalize(data$IL1A)
data$normalized_IRAK2 <- normalize(data$IRAK2)
data$normalized_PRKACB <- normalize(data$PRKACB)
data$normalized_TP53 <- normalize(data$TP53)

#for supressors

data$normalized_BIRC3 <- normalize1(data$BIRC3) 
data$normalized_TNFRSF10D <- normalize1(data$TNFRSF10D) 
data$normalized_PIK3CD <- normalize1(data$PIK3CD) 
data$normalized_PIK3R3 <- normalize1(data$PIK3R3) 
data$normalized_BCL2 <- normalize1(data$BCL2) 
data$normalized_BCL2L1 <- normalize1(data$BCL2L1) 
data$normalized_CAPN2 <- normalize1(data$CAPN2) 
data$normalized_AKT1 <- normalize1(data$AKT1) 


# adjusted score
data$Adjusted_Score <-  data$normalized_BCL2L11 + data$normalized_TP53INP1 + data$normalized_BIRC3 + data$normalized_CYCS + data$normalized_TNFRSF10D + data$normalized_PIK3CD + data$normalized_PIK3R3 + data$normalized_BCL2 + data$normalized_BCL2L1 + data$normalized_CAPN2 + data$normalized_PPP3CB + data$normalized_PRKX + data$normalized_ATM + data$normalized_IL1A + data$normalized_IRAK2 + data$normalized_PRKACB + data$normalized_TP53 + data$normalized_AKT1 # sum all them up

# category
data$category <- ifelse(data$Adjusted_Score > mean(data$Adjusted_Score), "more apoptotic", "less apoptotic")

#
print(data)


library('ggplot2')
ggplot(data, aes(x = group, y = Adjusted_Score, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("more apoptotic" = "red", "less apoptotic" = "blue")) +
  labs(x = "Group", y = "Adjusted Score", fill = "Category") +
  theme_minimal()


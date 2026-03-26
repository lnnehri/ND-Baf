#FA oxidation score
data <- data.frame(
  group = c(2, 3, 4),
  ABCD1 = c(1.103806896, NA, NA), # inducer
  ABCD3 = c(NA, NA, -1.1922394), # inducer
  ACAD11 = c(1.58162722, NA, 1.92399438), # inducer
  ACOT8 = c(NA, NA, 1.06839345), # inducer
  ACOX2 = c(1.514986236, 1.250884391, 1.661705853),
  AMACR = c(NA, NA, 1.11508322), # inducer
  CRAT = c(1.606008218, 1.215030617, 1.82216515), # inducer
  MCEE = c(NA, NA, 1.07960814), # inducer
  PCCB = c(-0.81918828, NA, -1.6613709), # inducer
  PEX13 = c(NA, NA, 1.06203812), # inducer
  PPARD = c(NA, NA, 1.09904979), # inducer
  SESN2 = c(1.93837683, 1.141918406, 2.255953786), # inducer
  SLC25A17 = c(-1.1810531, -1.036624976, -1.388890654) # inducer
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

data$normalized_ABCD1 <- normalize(data$ABCD1) # inducer doi: 10.1093/hmg/ddt645
data$normalized_ABCD3 <- normalize(data$ABCD3) # inducer https://doi.org/10.1093/hmg/ddu448
data$normalized_ACAD11 <- normalize(data$ACAD11) # inducer doi: 10.1016/j.ymgme.2010.12.005
data$normalized_ACOT8 <- normalize(data$ACOT8) # inducer https://doi.org/10.1016/j.bbadis.2012.03.009
data$normalized_ACOX2 <- normalize(data$ACOX2) # inducer https://doi.org/10.1073/pnas.1613228113
data$normalized_AMACR <- normalize(data$AMACR) # inducer DOI:10.1186/s40246-021-00344-1
data$normalized_CRAT <- normalize(data$CRAT) # inducer doi: 10.1016/j.molmet.2016.12.008
data$normalized_MCEE <- normalize(data$MCEE) # inducer doi: 10.1016/j.bbadis.2019.01.021
data$normalized_PCCB <- normalize(data$PCCB) # inducer doi: 10.3389/fgene.2022.807822
data$normalized_PEX13 <- normalize(data$PEX13) # inducer doi: 10.15252/embr.201642443
data$normalized_PPARD <- normalize(data$PPARD) # inducer DOI: 10.1016/j.bbrc.2009.12.127
data$normalized_SESN2 <- normalize(data$SESN2) # inducer doi: 10.1016/j.jcmgh.2021.04.015
data$normalized_SLC25A17 <- normalize(data$SLC25A17) # inducer doi: 10.3389/fcell.2020.00144

## bu genlerin inducer mi suppressor mu oldugunu bulman lazim
#hepsi inducer
data$Adjusted_Score <- data$normalized_ABCD1 + data$normalized_ABCD3 + data$normalized_ACAD11 + 
  data$normalized_ACOT8 + data$normalized_ACOX2 + data$normalized_AMACR + data$normalized_CRAT + 
  data$normalized_MCEE + data$normalized_PCCB + data$normalized_PEX13 + data$normalized_PPARD + 
  data$normalized_SESN2 + data$normalized_SLC25A17

data$category <- ifelse(data$Adjusted_Score > mean(data$Adjusted_Score), "more FA Oxidation", "less FA Oxidation")
print(data)

View(data)

for (i in 1:nrow(data)) {
  cat("Group:", data$group[i], "\t", 
      "Category:", data$category[i], "\n")
}

library('ggplot2')
ggplot(data, aes(x = group, y = Adjusted_Score, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("more FA Oxidation" = "red", "less FA Oxidation" = "blue")) +
  labs(x = "Group", y = "Adjusted Score", fill = "Category") +
  theme_minimal()


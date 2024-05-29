#ROS score
# df olusturma
data <- data.frame(
  group = c(2, 3, 4),
  AKT1 = c(NA, NA, -1.12537253),
  ATP7A = c(0.855987085, NA, 1.71298208),
  BCL2 = c(NA, NA, -1.4975039),
  BNIP3 = c(NA, 1.81056293, NA),
  CYBA = c(NA, NA, -1.9489035),
  CYP1A1 = c(NA, 1.08441409, NA),
  CYP1B1 = c(-1.5417554, NA, 2.46853471),
  DDIT4 = c(2.601605169, 2.961595164, 3.475183789),
  DHFR = c(-1.42470485, NA, -1.9538248),
  DUOX2 = c(1.985578193, 1.343795499, NA),
  EDN1 = c(-1.992393124, NA, -3.7231337),
  EGFR = c(1.07123324, NA, NA),
  EPHX2 = c(1.21480231, NA, 0.96828596),
  GCLM = c(-1.1789052, NA, NA),
  GPX1 = c(-1.080537207, NA, -1.1144963),
  GPX2 = c(NA, 1.37895444, NA),
  GPX3 = c(-1.150785991, NA, -1.2201409),
  HMOX1 = c(NA, NA, 1.20650602),
  HSP90AA1 = c(-1.945546836, -1.275173899, -2.087067788),
  MAOB = c(NA, NA, -2.7808074),
  NOS2 = c(-2.2962272, NA, -2.4195555),
  NOS3 = c(NA, 2.96993846, 3.33696695),
  NOSTRIN = c(-1.364725461, NA, -1.5761996),
  NOX1 = c(-2.98504751, -1.921708293, -4.35222649),
  NOXA1 = c(1.810805135, NA, NA),
  PARK7 = c(-0.924837713, NA, -1.1716096),
  PRDX1 = c(-1.469746004, NA, -1.8694652),
  PRDX3 = c(-0.7885859, NA, -1.1297303),
  PRDX4 = c(-0.977240804, NA, -1.3425915),
  PREX1 = c(1.322603536, 1.101901035, NA),
  PRKG2 = c(1.738295712, NA, 1.21980613),
  RFK = c(NA, NA, -1.4934019),
  RORA = c(1.174270541, NA, 1.21207093),
  SH3PXD2A = c(1.78368154, 1.88584993, 2.432630493),
  SH3PXD2B = c(NA, NA, -1.8178487),
  SLC7A2 = c(1.4209887, NA, 1.2661619),
  SOD1 = c(-1.123875745, NA, -1.497028),
  TXN = c(-1.455434486, NA, -1.8217603),
  TXNIP = c(-3.100828773, NA, -1.5381986),
  WASL = c(NA, NA, 1.00528093)
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


data$normalized_ATP7A <- normalize(data$ATP7A) # inducer https://doi.org/10.1038/s42003-020-0755-1
data$normalized_BNIP3 <- normalize(data$BNIP3) # inducer https://doi.org/10.1016/j.bbamcr.2015.02.022
data$normalized_CYBA <- normalize(data$CYBA) # inducer DOI: 10.1002/humu.21029
data$normalized_CYP1A1 <- normalize(data$CYP1A1) # inducer doi: 10.1016/j.cotox.2020.07.001
data$normalized_CYP1B1 <- normalize(data$CYP1B1) # inducer DOI: 10.1039/c9ra03674a
data$normalized_DDIT4 <- normalize(data$DDIT4) # inducer doi: 10.3389/fonc.2018.00106
data$normalized_DUOX2 <- normalize(data$DUOX2) # inducer doi: 10.1016/j.redox.2018.04.020
data$normalized_EDN1 <- normalize(data$EDN1) # inducer doi: 10.1038/sj.bjp.0706193
data$normalized_EGFR <- normalize(data$EGFR) # inducer doi: 10.1186/s13046-018-0728-0
data$normalized_EPHX2 <- normalize(data$EPHX2) # inducer https://doi.org/10.3389/fonc.2022.870721
data$normalized_HMOX1 <- normalize(data$HMOX1) # inducer DOI: 10.1089/ars.2005.7.80 
data$normalized_MAOB <- normalize(data$MAOB) # inducer doi: 10.1155/2017/3017947
data$normalized_NOS2 <- normalize(data$NOS2) # inducer https://doi.org/10.1016/j.niox.2019.04.009
data$normalized_NOX1 <- normalize(data$NOX1) # inducer doi: 10.1101/gad.339903.120
data$normalized_NOXA1 <- normalize(data$NOXA1) # inducer https://doi.org/10.1158/0008-5472.CAN-15-1512 
data$normalized_PREX1 <- normalize(data$PREX1) # inducer https://doi.org/10.4049/jimmunol.1002738
data$normalized_RFK <- normalize(data$RFK) # inducer https://doi.org/10.1016/j.yjmcc.2022.09.003
data$normalized_SH3PXD2A <- normalize(data$SH3PXD2A) # inducer https://doi.org/10.1242/jcs.203661
data$normalized_SH3PXD2B <- normalize(data$SH3PXD2B) # inducer https://doi.org/10.1242/jcs.203661
data$normalized_AKT1 <- normalize(data$AKT1) # inducer DOI:10.1126/scisignal.2003948
data$normalized_NOS3 <- normalize(data$NOS3) # inducer doi: 10.1172/JCI21968

data$normalized_WASL <- normalize1(data$WASL) # suppressor doi: 10.1016/j.parkreldis.2021.02.001
data$normalized_TXNIP <- normalize1(data$TXNIP) # suppressor https://doi.org/10.1038/s12276-023-01019-8
data$normalized_SLC7A2 <- normalize1(data$SLC7A2) # suppressor  DOI: https://doi.org/10.1124/dmd.121.000705
data$normalized_HSP90AA1 <- normalize1(data$HSP90AA1) # suppressor DOI:https://doi.org/10.1016/j.joca.2020.02.494
data$normalized_BCL2 <- normalize1(data$BCL2) # suppressor https://doi.org/10.1210/en.2015-1964
data$normalized_DHFR <- normalize1(data$DHFR) # suppressor doi: 10.1016/j.redox.2019.101185
data$normalized_GCLM <- normalize1(data$GCLM) # suppressor https://doi.org/10.1016/B978-0-12-420117-0.00001-3
data$normalized_GPX1 <- normalize1(data$GPX1) # suppressor doi: 10.3390/cancers14102560
data$normalized_GPX2 <- normalize1(data$GPX2) # suppressor DOI:10.1080/10715762.2021.1882677
data$normalized_GPX3 <- normalize1(data$GPX3) # suppressor https://doi.org/10.1016/j.bcp.2020.114365
data$normalized_NOSTRIN <- normalize1(data$NOSTRIN) # suppressor DOI:https://doi.org/10.1053/j.gastro.2006.12.035
data$normalized_PARK7 <- normalize1(data$PARK7) # suppressor https://doi.org/10.1016/j.exer.2019.107830
data$normalized_PRDX1 <- normalize1(data$PRDX1) # suppressor doi: 10.1186/s41021-021-00211-4
data$normalized_PRDX3 <- normalize1(data$PRDX3) # suppressor DOI:10.3892/etm.2021.9870
data$normalized_PRDX4 <- normalize1(data$PRDX4) # suppressor doi: 10.1177/1533033819864313
data$normalized_PRKG2 <- normalize1(data$PRKG2) # suppressor doi: 10.1016/j.ajpath.2016.10.016
data$normalized_RORA <- normalize1(data$RORA) # suppressor  https://doi.org/10.3390/ijms221910665
data$normalized_SOD1 <- normalize1(data$SOD1) # suppressor doi: 10.1016/j.jneuroim.2006.10.003
data$normalized_TXN <- normalize1(data$TXN) # suppressor DOI: 10.1080/10245332.2016.1173341

###

data$Adjusted_Score <- data$normalized_ATP7A + data$normalized_BNIP3 + data$normalized_CYBA + 
  data$normalized_CYP1A1 + data$normalized_CYP1B1 + data$normalized_DDIT4 + 
  data$normalized_DUOX2 + data$normalized_EDN1 + data$normalized_EGFR + 
  data$normalized_EPHX2 + data$normalized_HMOX1 + data$normalized_HSP90AA1 + 
  data$normalized_MAOB + data$normalized_NOS2 + data$normalized_NOX1 + 
  data$normalized_NOXA1 + data$normalized_PREX1 + data$normalized_RFK + 
  data$normalized_SH3PXD2A + data$normalized_SH3PXD2B + data$normalized_SLC7A2 + 
  data$normalized_TXNIP + data$normalized_WASL + data$normalized_AKT1 + 
  data$normalized_BCL2 + data$normalized_DHFR + data$normalized_GCLM + 
  data$normalized_GPX1 + data$normalized_GPX2 + data$normalized_GPX3 + 
  data$normalized_NOS3 + data$normalized_NOSTRIN + data$normalized_PARK7 + 
  data$normalized_PRDX1 + data$normalized_PRDX3 + data$normalized_PRDX4 + 
  data$normalized_PRKG2 + data$normalized_RORA + data$normalized_SOD1 + 
  data$normalized_TXN
# ROS ayrimi
data$category <- ifelse(data$Adjusted_Score > mean(data$Adjusted_Score), "more ROS", "less ROS")
print(data)

View(data)

for (i in 1:nrow(data)) {
  cat("Group:", data$group[i], "\t", 
      "Category:", data$category[i], "\n")
}

library('ggplot2')
ggplot(data, aes(x = group, y = Adjusted_Score, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("more ROS" = "red", "less ROS" = "blue")) +
  labs(x = "Group", y = "Adjusted Score", fill = "Category") +
  theme_minimal()

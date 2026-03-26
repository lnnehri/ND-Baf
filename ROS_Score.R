#ROS score
#synergy R code 
#for p_adj<0.01 and unique ROS genes 
library(readr)
library(dplyr)
data_ROS <- read_delim("Desktop/data_ROS.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(data_ROS)

#adding character column
data_ROS$char<-NA
View(data_ROS)

###ROS genes
data_ROS <- data_ROS %>%
  mutate(char = case_when(
    is.na(`ND+BAF`) ~ "synergy",
    is.na(ND) & is.na(Baf) ~ "synergy",
    is.na(Baf) & (ND > 0) & (ND + 0.4 < `ND+BAF`) ~ "synergy",
    is.na(Baf) & (ND > 0) & (ND > `ND+BAF` + 0.4) ~ "synergy",
    is.na(Baf) & (ND < 0) & (ND - 0.3 < `ND+BAF`) ~ "synergy",
    is.na(Baf) & (ND < 0) & (ND > `ND+BAF` - 0.3) ~ "synergy",
    is.na(ND) & (Baf > 0) & (Baf + 0.4 < `ND+BAF`) ~ "synergy",
    is.na(ND) & (Baf > 0) & (Baf > `ND+BAF` - 0.3) ~ "synergy",
    is.na(ND) & (Baf < 0) & (Baf - 0.3 < `ND+BAF`) ~ "synergy",
    is.na(ND) & (Baf < 0) & (Baf > `ND+BAF` - 0.3) ~ "synergy",
    (ND > 0) & (Baf > 0) & ((ND * Baf) + 0.4 < `ND+BAF`) ~ "synergy",
    (ND > 0) & (Baf > 0) & ((ND * Baf) > (`ND+BAF` + 0.4)) ~ "synergy",
    (ND < 0) & (Baf < 0) & (-(ND * Baf) - 0.3 < `ND+BAF`) ~ "synergy",
    (ND < 0) & (Baf < 0) & (-(ND * Baf) > (`ND+BAF` - 0.3)) ~ "synergy",
    TRUE ~ "additive"
  ))

View(data_ROS)

#adding regulation column
data_ROS$regul<-NA
View(data_ROS)

#filling increase or decrease characteristics of genes
data_ROS <- data_ROS %>%
  mutate(regul = case_when(
    char == "additive" & (`ND+BAF` > ND) & is.na(Baf) ~ "increase",
    char == "additive" & (`ND+BAF` < ND) & is.na(Baf) ~ "decrease",
    char == "additive" & (`ND+BAF` < Baf ) & is.na(ND)~ "decrease",
    char == "additive" & (`ND+BAF` > Baf ) & is.na(ND)~ "increase",
    is.na(ND) & is.na(Baf) & `ND+BAF` > 0 ~ "increase",
    is.na(ND) & is.na(Baf) & `ND+BAF` < 0 ~ "decrease",
    is.na(Baf) & ND > 0 & ND < `ND+BAF` ~ "increase",
    is.na(Baf) & ND > 0 & ND > `ND+BAF` ~ "decrease",
    is.na(Baf) & ND < 0 & ND < `ND+BAF` ~ "increase",
    is.na(Baf) & ND < 0 & ND > `ND+BAF` ~ "decrease",
    is.na(ND) & Baf > 0 & Baf < `ND+BAF` ~ "increase",
    is.na(ND) & Baf > 0 & Baf > `ND+BAF` ~ "decrease",
    is.na(ND) & Baf < 0 & Baf < `ND+BAF` ~ "increase",
    is.na(ND) & Baf < 0 & Baf > `ND+BAF` ~ "decrease",
    is.na(ND) & Baf > 0 & is.na(`ND+BAF`) ~ "decrease",
    is.na(ND) & Baf < 0 & is.na(`ND+BAF`) ~ "increase",
    is.na(Baf) & ND < 0 & is.na(`ND+BAF`) ~ "increase",
    is.na(Baf) & ND > 0 & is.na(`ND+BAF`) ~ "decrease",
    is.na(`ND+BAF`) & ND > 0 & Baf > 0 ~ "decrease",
    is.na(`ND+BAF`) & ND < 0 & Baf < 0 ~ "increase",
    ND > 0 & Baf > 0 & (ND * Baf) + 0.4 < `ND+BAF` ~ "increase",
    ND > 0 & Baf > 0 & (ND * Baf) > `ND+BAF` + 0.4 ~ "decrease",
    char == "additive" & ND > 0 & Baf > 0 & (ND * Baf) < `ND+BAF` ~ "increase",
    char == "additive" & ND > 0 & Baf > 0 & (ND * Baf) > `ND+BAF` ~ "decrease",
    ND < 0 & Baf < 0 & -(ND * Baf) - 0.3 < `ND+BAF` ~ "increase",
    ND < 0 & Baf < 0 & -(ND * Baf) > `ND+BAF` - 0.3 ~ "decrease",
    char == "additive" &  ND < 0 & Baf < 0 & -(ND * Baf) < `ND+BAF` ~ "increase",
    char == "additive" &  ND < 0 & Baf < 0 & -(ND * Baf) > `ND+BAF`~ "decrease",
    TRUE ~ NA_character_
  ))

View(data_ROS)

#adding gene up or down column column
data_ROS$up_down<-NA

data_ROS <- data_ROS %>%
  mutate(up_down = case_when(
    ND > 0 | Baf > 0 | `ND+BAF` > 0 ~ "up",
    TRUE ~ "down"
  ))


View(data_ROS)


#exploring desired genes
unique(data_ROS$regul)
unique(data_ROS$up_down)
unique(data_ROS$char)

data_ROS$symbol[data_ROS$regul == "increase" & data_ROS$up_down == "down" & data_ROS$char == "synergy"]

print(data_ROS)

#Code part is finished
#You can see the table below
GENERATED ROS TABLE:

entrezgene	Gene	ND	Baf	ND+BAF	char	regul	up_down
207	AKT1	NA	NA	-1.12537253	synergy	decrease	down
538	ATP7A	0.855987085	NA	1.71298208	synergy	increase	up
596	BCL2	NA	NA	-1.4975039	synergy	decrease	down
664	BNIP3	NA	1.81056293	NA	synergy	decrease	up
1535	CYBA	NA	NA	-1.9489035	synergy	decrease	down
1543	CYP1A1	NA	1.08441409	NA	synergy	decrease	up
1545	CYP1B1	-1.5417554	NA	2.46853471	synergy	increase	up
54541	DDIT4	2.601605169	2.961595164	3.475183789	synergy	decrease	up
1719	DHFR	-1.42470485	NA	-1.9538248	synergy	decrease	down
50506	DUOX2	1.985578193	1.343795499	NA	synergy	decrease	up
1906	EDN1	-1.992393124	NA	-3.7231337	synergy	decrease	down
1956	EGFR	1.07123324	NA	NA	synergy	decrease	up
2053	EPHX2	1.21480231	NA	0.96828596	additive	decrease	up
2730	GCLM	-1.1789052	NA	NA	synergy	increase	down
2876	GPX1	-1.080537207	NA	-1.1144963	synergy	decrease	down
2877	GPX2	NA	1.37895444	NA	synergy	decrease	up
2878	GPX3	-1.150785991	NA	-1.2201409	synergy	decrease	down
3162	HMOX1	NA	NA	1.20650602	synergy	increase	up
3320	HSP90AA1	-1.945546836	-1.275173899	-2.087067788	synergy	increase	down
4129	MAOB	NA	NA	-2.7808074	synergy	decrease	down
4843	NOS2	-2.2962272	NA	-2.4195555	synergy	decrease	down
4846	NOS3	NA	2.96993846	3.33696695	additive	increase	up
115677	NOSTRIN	-1.364725461	NA	-1.5761996	synergy	decrease	down
27035	NOX1	-2.98504751	-1.921708293	-4.35222649	synergy	increase	down
10811	NOXA1	1.810805135	NA	NA	synergy	decrease	up
11315	PARK7	-0.924837713	NA	-1.1716096	synergy	decrease	down
5052	PRDX1	-1.469746004	NA	-1.8694652	synergy	decrease	down
10935	PRDX3	-0.7885859	NA	-1.1297303	synergy	decrease	down
10549	PRDX4	-0.977240804	NA	-1.3425915	synergy	decrease	down
57580	PREX1	1.322603536	1.101901035	NA	synergy	decrease	up
5593	PRKG2	1.738295712	NA	1.21980613	synergy	decrease	up
55312	RFK	NA	NA	-1.4934019	synergy	decrease	down
6095	RORA	1.174270541	NA	1.21207093	additive	increase	up
9644	SH3PXD2A	1.78368154	1.88584993	2.432630493	synergy	decrease	up
285590	SH3PXD2B	NA	NA	-1.8178487	synergy	decrease	down
6542	SLC7A2	1.4209887	NA	1.2661619	additive	decrease	up
6647	SOD1	-1.123875745	NA	-1.497028	synergy	decrease	down
7295	TXN	-1.455434486	NA	-1.8217603	synergy	decrease	down
10628	TXNIP	-3.100828773	NA	-1.5381986	synergy	increase	down
8976	WASL	NA	NA	1.00528093	synergy	increase	up

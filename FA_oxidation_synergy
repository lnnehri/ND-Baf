#FA_oxidation score
#synergy R code 
#for p_adj<0.01 and unique FA_oxidation
library(readr)
library(dplyr)
data_FA_oxidation <- read_delim("Desktop/FA_oxidation.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(data_FA_oxidation)

#adding character column
data_FA_oxidation$char<-NA
View(data_FA_oxidation)

###FA_oxidation genes
data_FA_oxidation <- data_FA_oxidation %>%
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

View(data_FA_oxidation)

#adding regulation column
data_FA_oxidation$regul<-NA
View(data_FA_oxidation)

#filling increase or decrease characteristics of genes
data_FA_oxidation <- data_FA_oxidation %>%
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

View(data_FA_oxidation)

#adding gene up or down column column
data_FA_oxidation$up_down<-NA

data_FA_oxidation <- data_FA_oxidation %>%
  mutate(up_down = case_when(
    ND > 0 | Baf > 0 | `ND+BAF` > 0 ~ "up",
    TRUE ~ "down"
  ))


View(data_FA_oxidation)


#exploring desired genes
unique(data_FA_oxidation$regul)
unique(data_FA_oxidation$up_down)
unique(data_FA_oxidation$char)

data_FA_oxidation$symbol[data_FA_oxidation$regul == "increase" & data_FA_oxidation$up_down == "down" & data_FA_oxidation$char == "synergy"]

print(data_FA_oxidation)

##codei is finished
##the table below:
GENE ID	GENE	ND	Baf	ND+BAF	char	regul	up_down
215	ABCD1	1.103806896	NA	NA	synergy	decrease	up
5825	ABCD3	NA	NA	-1.1922394	synergy	decrease	down
84129	ACAD11	1.58162722	NA	1.92399438	additive	increase	up
10005	ACOT8	NA	NA	1.06839345	synergy	increase	up
8309	ACOX2	1.514986236	1.250884391	1.661705853	additive	decrease	up
23600	AMACR	NA	NA	1.11508322	synergy	increase	up
1384	CRAT	1.606008218	1.215030617	1.82216515	additive	decrease	up
84693	MCEE	NA	NA	1.07960814	synergy	increase	up
5096	PCCB	-0.81918828	NA	-1.6613709	synergy	decrease	down
5194	PEX13	NA	NA	1.06203812	synergy	increase	up
5467	PPARD	NA	NA	1.09904979	synergy	increase	up
83667	SESN2	1.93837683	1.141918406	2.255953786	additive	increase	up
10478	SLC25A17	-1.1810531	-1.036624976	-1.388890654	synergy	increase	down

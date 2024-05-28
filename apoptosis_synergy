#apoptosis synergy R code 
#for p_adj<0.01 and apoptosis genes from table xx. of manuscript
library(readr)
library(dplyr)
data_apop <- read_delim("Desktop/data_apop.csv", delim = ";", 
                   +     escape_double = FALSE, trim_ws = TRUE)
View(data_apop)

#adding character column
data_apop$char<-NA
View(data_apop)

###apoptosis genes
data_apop <- data_apop %>%
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

View(data_apop)

#adding regulation column
data_apop$regul<-NA
View(data_apop)

#filling increase or decrease characteristics of genes
data_apop <- data_apop %>%
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

View(data_apop)

#adding gene up or down column column
data_apop$up_down<-NA

data_apop <- data_apop %>%
  mutate(up_down = case_when(
    ND > 0 | Baf > 0 | `ND+BAF` > 0 ~ "up",
    TRUE ~ "down"
  ))


View(data_apop)


#exploring desired genes
unique(data_apop$regul)
unique(data_apop$up_down)
unique(data_apop$char)

data_apop$symbol[data_apop$regul == "increase" & data_apop$up_down == "down" & data_apop$char == "synergy"]

print(data_apop)
# A tibble: 19 × 8
#id symbol       ND    Baf `ND+BAF` char     regul    up_down
#<dbl> <chr>     <dbl>  <dbl>    <dbl> <chr>    <chr>    <chr>  
#  1 10018 BCL2L11    1.64  0.945     1.71 additive increase up     
#2 94241 TP53INP1   1.61  1.42      1.75 synergy  decrease up     
#3  5293 PIK3CD    NA     1.41     NA    synergy  decrease up     
#4  8503 PIK3R3     1.95 NA         1.91 additive decrease up     
#5   472 ATM       NA    NA         1.02 synergy  increase up     
#6  3552 IL1A      NA    NA         1.84 synergy  increase up     
#7  3656 IRAK2     NA    NA         1.11 synergy  increase up     
#8  5567 PRKACB    NA    NA         1.28 synergy  increase up     
#9  7157 TP53      NA    NA         1.18 synergy  increase up     
#10   596 BCL2      NA    NA        -1.50 synergy  decrease down   
#11   598 BCL2L1    NA    NA        -1.62 synergy  decrease down   
#12   824 CAPN2     NA    NA        -1.19 synergy  decrease down   
#13  5532 PPP3CB    NA    NA        -1.60 synergy  decrease down   
#14  5613 PRKX      NA    NA        -1.74 synergy  decrease down   
#15   207 AKT1      NA    NA        -1.13 synergy  decrease down   
#16   330 BIRC3     -1.29 NA        -1.08 synergy  increase down   
#17 54205 CYCS      -1.54 NA        NA    synergy  increase down   
#18  8793 TNFRSF10D -1.05 NA        NA    synergy  increase down   
#19 23052 ENDOD1    -1.36 -0.852    -1.89 synergy  decrease down  

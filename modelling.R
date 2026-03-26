################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

##USE 
#SIGMOID FOR DEATH PATTERN OF CELLS

# sigmoid function creation
sigmoid <- function(x, L, x0, k) {
  L / (1 + exp(k * (x - x0)))
}

x <- seq(0, 96, length.out = 100)
# sigmoid parameters
L <- 45000  # max number of cells
x0 <- 72    # middle point for acceleration of death
k <- 0.1    # for the curve

# use sigmoid function to calculate y values
y <- sigmoid(x, L, x0, k)

data <- data.frame(x = x, y = y)

ggplot(data, aes(x, y)) +
  geom_line(color = "blue", size = 1.15) +
  geom_vline(xintercept = c(48, 72, 96), linetype = "dashed", color = "red") +
  geom_text(data = data.frame(x = c(48, 72,96), y = 45000, label = c("48. hour", "72. hour", '96.hour')), aes(label = label), vjust = 0, hjust = 1.1, color = "red", fontface = "bold") +
  labs(title = "Death Pattern",
       x = "Time (hour)",
       y = "Number of Living Cells") +
  annotate("text", x = 24, y = 25000, label = "Metastatic Phenotype Change", color = "red", size = 4, fontface = "bold") +
  annotate("rect", xmin = 6, xmax = 42, ymin = 23850, ymax = 26000, fill = "transparent", color = "red") +
  
  
  theme_minimal() +
  theme(
    text = element_text(face = "bold"), # bold
    axis.title.y = element_text(margin = margin(r = 20)) # for y
  )

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################


### USE
#POISSON OF METASTATIC/DEATH/ROUNDED:

#shape change probabilities:
lambda <- 26 / 126 #for metastatic cells in all living cells
lambda
#prob of seeing a death cell until 48 hours:
prob_death<-dpois(1,0.12) #lamda is calculated from 12 cells are death among 100 cells,88 cells living
#seeing a living cell until 48 hours:
prob_live<-1-prob_death
#among living cells, phenotypic change:
prob_metastasis<-dpois(1, lambda) 
prob_rounded<- 1-prob_metastasis
#living and metastatic and living and rounded:
living_metastatic<-prob_metastasis*prob_live
living_rounded<-prob_rounded*prob_live
#cheking total prob is 1 or not:
living_metastatic + living_rounded + prob_death

#result is 1, then plotting it:
#data frame creation
data <- data.frame(Shape = c("Rounded", "Metastatic", "Death"),
                   Probability = c(0.7435614, 0.1500081, 0.1064305))

#pie chart
ggplot(data, aes(x = "", y = Probability, fill = Shape)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  labs(title = "Phenotypic Distribution of Cells Until 48. Hour",
       x = NULL,
       y = NULL,
  ) +
  scale_fill_manual(values = c("Rounded" = "blue", "Metastatic" = "red", "Death" = "green")) +
  geom_text(aes(label = paste(round(Probability * 100), "%")), position = position_stack(vjust = 0.5), color = "black", size = 8) +
  theme_void() +
  theme( legend.position = "right",
         plot.title = element_text(face = "bold", size = 15, margin = margin(b = 1)), 
         axis.title = element_text(face = "bold"), # lables
         legend.text = element_text(face = "bold"), # lejant
         plot.margin = margin(t = 10, b = 20) 
  )

#THIS POISSON FINISHED


################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

##BU ALTTAKINDE 3 GAMMA DISTRIBUTIONU ALT ALTA YAPIYORSUN
##VERSION 1
#BUNA DONUP BAK 
#asagida normalization yapcan
###########
library(cowplot)
#asagidaki 100-500den bakinca da farkli cikiyor
#birinden birisi sabit olacak
#literaturle destekle altta distributionlari
#oleni ele, FA cok olani sonra 
# x de??erleri
x <- seq(0, 100, length.out = 500)
#lamda values -- normalization yapilacak sonra ayni seye getireceksin 
#scaleler ayni olmasi lazim

#bir katsayiyla ortak bir formuluzasyon etmen lazim

# Lambda ve a de??erleri
lambda_values <- c(prob_death, prob_live, prob_live) #lipid ratios
a_values <- c(prob_death, prob_rounded, prob_metastasis) #calcium

# Veri ??er??evesini olu??turma
data <- expand.grid(x = x, lambda = lambda_values, a = a_values)
data$y <- dgamma(data$x, shape = data$a, rate = data$lambda)

# Grafiklerin olu??turulmas??
plot1 <- ggplot(data[data$lambda == prob_death,], aes(x = x, y = y)) +
  geom_line(color = "green") +
  labs(title = "Death (Lambda=prob_death)", x = "FA concentration", y = "Ca+ concentration") +
  theme_minimal()

plot2 <- ggplot(data[data$lambda == prob_live,], aes(x = x, y = y)) +
  geom_line(color = "blue") +
  labs(title = "Rounded (Lambda=prob_live)", x = "FA concentration", y = "Ca+ concentration") +
  theme_minimal()

plot3 <- ggplot(data[data$lambda == prob_live,], aes(x = x, y = y)) +
  geom_line(color = "red") +
  labs(title = "Metastatic (Lambda=prob_live)", x = "FA concentration", y = "Ca+ concentration") +
  theme_minimal()

# ???? grafi??i ??st ??ste birle??tirme
final_plot <- plot_grid(plot1, plot2, plot3, nrow = 3)

final_plot
##yukaridaki bitti
#yukaridaki elde edilen grafikte, rounded ve metastatik arasinda anlamli bir fark var mi yok mu bakmak icin
#t test yapiyoruz:
y <- rgamma(data$x, shape = prob_death, rate = prob_death) #dead #green
z <- rgamma(data$x, shape = prob_live, rate = prob_rounded) #rounded #blue
x <- rgamma(data$x, shape = prob_live, rate = prob_metastasis) #metastasis #red
t.test(x,z)

###literaturden destek lazim bos sallamamak icin
#yukaridakini yessili elemek icin kyllancan
#bunu metastatigi almak icin

data <- data.frame(
  group = c(rep("x", length(x)), rep("z", length(z)), rep("y", length(y))),
  value = c(x, z, y)
)
anova_result <- aov(value ~ group, data = data)
print(summary(anova_result))
ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  labs(title = "Gruplar Arasındaki Ortalama Farkın Karşılaştırılması",
       x = "Grup",
       y = "Değer") +
  theme_minimal()
#####


################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

#bu opsiyonda direkt farkli bir sekilde metastatikin gammasini modelliyorsun:

######su da 2. opsiyon olarak donup bak:
###odagin bu su an :::::
###probability ofDEATH OR LIVING BASED ON LIPID ACCUMULATION

L_treshold<- prob_death/(prob_death+prob_live) #bunu wetlabdan gelene gore guncelle - bu deger fix olacak
L_treshold<- prob_death
L_hour<- prob_live/(prob_death+prob_live)
L_hour

#bunu wetlabdan gelene gore guncelle - bu deger degisken olabilir
#The pgamma function is the probability of waiting the 48 hours until we have the arrival of alpha=L_treshold of lipids in a cell. 
#Here, we assume that L_hour of lipids accumulate in the cell every hour. 
#Here we calculate the probability that we have to wait 48 hours to get L_treshold of lipid accumulation.
#we are using x=0.5, to represent the halftime of the cells (48/96)

lipid_accumulation_prob<-dgamma(0.5, shape=(L_treshold),scale=1/(L_hour))
lipid_accumulation_prob

#deneme
# Gama dağılımı parametreleri
shape <- L_treshold  # Şekil parametresi
rate <- lipid_accumulation_prob # Oran parametresi

# Olasılık yoğunluk fonksiyonunu hesapla
x <- seq(0, 100, length=100) # 0'dan 100'a kadar x değerlerini oluştur
y <- dgamma(x, shape, rate) # Gama dağılımının olasılık yoğunluk fonksiyonunu hesapla

# Grafiği çiz
plot(x, y, type="l", col="blue", lwd=2, 
     main="Gamma Dağılımı", xlab="x", ylab="Olasılık Yoğunluğu")

#ustteki 2. opsiyon bitti

###odagin bu su an :::::
###ROUNDED OR METASTATIC BASED ON CALCIUM ACCUMULATION
#Here we calculate the probability that we have to wait 48 hours to get Ca_treshold of calcium accumulation.
#The pgamma function is the probability of waiting the 48 hours until we have the arrival of alpha=ca_treshold of calciums in a living cell. 
#Here, we assume that ca_hour of lipids accumulate in the living cell every hour. 

Ca_treshold<- prob_metastasis /(prob_live) #update et - bu deger fix olacak
Ca_hour<- prob_live #update et - bu deger degisken olabilir
calcium_accumulation_prob<-dgamma(0.5, shape=Ca_treshold,scale=1/Ca_hour)
calcium_accumulation_prob

#deneme
# Gama dağılımı parametreleri
shape <- calcium_accumulation_prob  # Şekil parametresi
rate <- lipid_accumulation_prob # Oran parametresi

# Olasılık yoğunluk fonksiyonunu hesapla
x <- seq(0, 100, length=100) # 0'dan 100'a kadar x değerlerini oluştur
y <- dgamma(x, shape, rate) # Gama dağılımının olasılık yoğunluk fonksiyonunu hesapla

# Grafiği çiz
plot(x, y, type="l", col="blue", lwd=2, 
     main="Gamma Dağılımı", xlab="x", ylab="Olasılık Yoğunluğu")

#ustteki 2. opsiyon bitti



################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

#bu 3.de de ortalama kac saatte ulasacak bu metastatik karar alma mekanizmasina ona bakiyorsun:

#odagin bu su. an::::
#########DECIDING METASTATIC OR NOT BASED ON LIPID AND CALCIUM
#shape=calcium level, scale=1/fatty acid level
#buradaki L_current belli bir tresholdun uzerindeki 
#bu degerleri bir birine gore normalize etmen gerekiyor oncesinde
#bastaki kumulatif dagilim degerine ne verecegin hala cok guzel ve buyuk bir soru isareti
#ona normalize ettikten sonra degerleri bakabilirsin belki 
#pgamma belli bir degerin altindaki kumulatif olasiligi hesaplamak icin kullaniliyor
#dgamma belli bir x degerinin gelme olasiligini hesaplar
#ben burada x degerine odaklanabilirim 
#burada soyle bir kod yazabilirsin:
#bir hucren var, onun L-hour ve Ca_hour valuelerini bildigini varsayiyorum:
#bu hucre icin sunlari hesapliyorsun:
#if lipid accumulation probability is greater than x yani bir prob tayin ediyorsun
#sonra and if calcium accumulation probabilty is greater than y 
# eger bu kosullara uyuyorsa, kacinci saatte bu degerlere ulasacagini hem ca hem FA icin hesapla

#deneme
#Ca_current<-2000
#L_current<-40
#dgamma(600, shape=Ca_current,scale=1/L_current)

#3. odagin bu su an:
##bu alttaki guzel::
#ya da direkt chatgtpden aldim suna bakarsin modifye edip:
find_threshold_time <- function(accumulation_value, threshold_value, shape, scale) {
  time <- 0
  cumulative_probability <- 0
  
  while (cumulative_probability < threshold_value) {
    time <- time + 1
    cumulative_probability <- pgamma(accumulation_value * time, shape = shape, scale = scale)
  }
  
  return(time)
}
#shape=calcium level, scale=1/fatty acid level
# Kullanım örneği:
# Saatlik birikim değeri
#treshold value icin deger belirlemek lazim 
#accumulation value icin deger belirlremek lazim
accumulation_value <- calcium_accumulation_prob
# Treshold değeri 100 birim birikim olsun
threshold_value <- prob_live
# Gamma dağılımı parametreleri
shape <- calcium_accumulation_prob
scale <- 1/lipid_accumulation_prob

# Treshold değerine ulaşma süresini bulma
threshold_time <- find_threshold_time(accumulation_value, threshold_value, shape, scale)
print(paste("Treshold değerine", threshold_value, "ulaşma süresi:", threshold_time, "saat"))

###grafigi

# Birikim değeri
time_values <- seq(0, 100, by = 1)
cumulative_prob_values <- pgamma(accumulation_value * time_values, shape = shape, scale = scale)

# Veri çerçevesi oluşturma
df <- data.frame(Time = time_values, Cumulative_Probability = cumulative_prob_values)

# Grafik oluşturma
ggplot(df, aes(x = Time, y = Cumulative_Probability)) +
  geom_line() +
  geom_hline(yintercept = threshold_value, linetype = "dashed", color = "red") +
  geom_vline(xintercept = threshold_time, linetype = "dashed", color = "blue") +
  labs(x = "Time", y = "Cumulative Probability", title = "Cumulative Probability Over Time") +
  theme_minimal()

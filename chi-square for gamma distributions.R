# Gerekli paketleri yükleyelim
if(!require("ggplot2")) install.packages("ggplot2", dependencies=TRUE)
library(ggplot2)

# Parametreler
prob_death <- 0.1064305
prob_live <- 0.8935695
prob_rounded <- 0.7435614
prob_metastasis <- 0.1500081

# Gamma dağılımlarından veriler oluşturma
set.seed(123)
x <- rgamma(500, shape = prob_live, rate = prob_metastasis) # Metastatik
z <- rgamma(500, shape = prob_live, rate = prob_rounded)    # Yuvarlak
y <- rgamma(500, shape = prob_death, rate = prob_death)     # Ölüm

# Sürekli verileri kategorik hale getirmek için binleme işlemi
breaks <- seq(0, 50, by = 10) # 5 kategori
labels <- c("0-10", "10-20", "20-30", "30-40", "40-50")

x_cat <- cut(x, breaks = breaks, labels = labels, include.lowest = TRUE)
z_cat <- cut(z, breaks = breaks, labels = labels, include.lowest = TRUE)
y_cat <- cut(y, breaks = breaks, labels = labels, include.lowest = TRUE)

# Gözlenen frekansları hesapla
obs_x <- table(x_cat)
obs_z <- table(z_cat)
obs_y <- table(y_cat)

# Beklenen frekansları hesapla (eşit dağılım varsayıyoruz)
observed <- rbind(obs_x, obs_z, obs_y)
expected <- colMeans(observed)

# Chi-Square testi
chisq_x <- chisq.test(obs_x, p = expected / sum(expected))
chisq_z <- chisq.test(obs_z, p = expected / sum(expected))
chisq_y <- chisq.test(obs_y, p = expected / sum(expected))

# Sonuçları gösterelim
cat("Metastatik Grup için Chi-Square Sonuçları:\n")
print(chisq_x)

cat("\nYuvarlak Grup için Chi-Square Sonuçları:\n")
print(chisq_z)

cat("\nÖlüm Grubu için Chi-Square Sonuçları:\n")
print(chisq_y)

# Grafik için veriyi hazırlayalım
data <- data.frame(
  group = c(rep("x", length(x_cat)), rep("z", length(z_cat)), rep("y", length(y_cat))),
  category = c(as.character(x_cat), as.character(z_cat), as.character(y_cat))
)

# Grupların kategorik dağılımını gösteren barplot
ggplot(data, aes(x = category, fill = group)) +
  geom_bar(position = "dodge") +
  labs(title = "Grupların Kategorik Dağılımı",
       x = "Kategori",
       y = "Frekans") +
  theme_minimal() +
  theme(
    text = element_text(face = "bold"), # Kalın yazılar
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 15)
  )
##########################

# Load necessary packages
if(!require("ggplot2")) install.packages("ggplot2", dependencies=TRUE)
library(ggplot2)

# Parameters for gamma distributions
prob_death <- 0.1064305
prob_live <- 0.8935695
prob_rounded <- 0.7435614
prob_metastasis <- 0.1500081

# Generate data from gamma distributions
set.seed(123)
x <- rgamma(500, shape = prob_live, rate = prob_metastasis) # Metastatic
z <- rgamma(500, shape = prob_live, rate = prob_rounded)    # Rounded
y <- rgamma(500, shape = prob_death, rate = prob_death)     # Death

# Categorize continuous data into bins
breaks <- seq(0, 50, by = 10) # 5 categories
labels <- c("0-10", "10-20", "20-30", "30-40", "40-50")

x_cat <- cut(x, breaks = breaks, labels = labels, include.lowest = TRUE)
z_cat <- cut(z, breaks = breaks, labels = labels, include.lowest = TRUE)
y_cat <- cut(y, breaks = breaks, labels = labels, include.lowest = TRUE)

# Calculate observed frequencies
obs_x <- table(x_cat)
obs_z <- table(z_cat)
obs_y <- table(y_cat)

# Calculate expected frequencies (assuming equal distribution)
observed <- rbind(obs_x, obs_z, obs_y)
expected <- colMeans(observed)

# Chi-Square tests
chisq_x <- chisq.test(obs_x, p = expected / sum(expected))
chisq_z <- chisq.test(obs_z, p = expected / sum(expected))
chisq_y <- chisq.test(obs_y, p = expected / sum(expected))

# Display results
cat("Chi-Square Results for Metastatic Group:\n")
print(chisq_x)

cat("\nChi-Square Results for Rounded Group:\n")
print(chisq_z)

cat("\nChi-Square Results for Death Group:\n")
print(chisq_y)

# Prepare data for visualization
data <- data.frame(
  group = c(rep("Metastatic", length(x_cat)), rep("Rounded", length(z_cat)), rep("Death", length(y_cat))),
  category = c(as.character(x_cat), as.character(z_cat), as.character(y_cat))
)

# Barplot showing categorical distributions of groups
ggplot(data, aes(x = category, fill = group)) +
  geom_bar(position = "dodge") +
  labs(title = "Categorical Distribution of Groups",
       x = "Category",
       y = "Frequency") +
  theme_minimal() +
  theme(
    text = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 15)
  )
###########################################################

# Calculate residuals
residuals <- observed - matrix(rep(expected, nrow(observed)), nrow = nrow(observed), byrow = TRUE)

# Convert to data frame for visualization
res_data <- data.frame(
  Group = rep(c("Metastatic", "Rounded", "Death"), each = length(labels)),
  Category = rep(labels, times = 3),
  Residual = as.vector(residuals)
)

# Plotting residuals
ggplot(res_data, aes(x = Category, y = Residual, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Residuals from Chi-Square Test",
       x = "Category",
       y = "Residual (Observed - Expected)") +
  theme_minimal() +
  theme(
    text = element_text(face = "bold"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 15)
  )

############################

# Gerekli paketleri yükleyelim
if(!require("vcd")) install.packages("vcd", dependencies=TRUE)
library(vcd)

# Gözlenen frekansları birleştirelim
observed <- rbind(obs_x, obs_z, obs_y)
colnames(observed) <- labels
rownames(observed) <- c("Metastatic", "Rounded", "Death")

# Mosaic Plot
mosaic(observed, 
       shade = TRUE, 
       legend = TRUE,
       main = "Mosaic Plot of Observed Frequencies",
       labeling = labeling_values)

########################################

# Residuals hesapla
residuals <- observed - matrix(rep(expected, nrow(observed)), nrow = nrow(observed), byrow = TRUE)

# Residuals verisini düzenleyelim
res_data <- data.frame(
  Group = rep(c("Metastatic", "Rounded", "Death"), each = length(labels)),
  Category = rep(labels, times = 3),
  Residual = as.vector(residuals)
)

# Heatmap çizimi
ggplot(res_data, aes(x = Category, y = Group, fill = Residual)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(min(res_data$Residual), max(res_data$Residual)),
                       space = "Lab", name="Residuals") +
  labs(title = "Heatmap of Chi-Square Residuals",
       x = "Category",
       y = "Group") +
  theme_minimal()

############################


# Gözlenen ve beklenen frekansları aynı tabloya koy
expected_data <- data.frame(
  Group = rep(c("Metastatic", "Rounded", "Death"), each = length(labels)),
  Category = rep(labels, times = 3),
  Frequency = c(expected, expected, expected),
  Type = rep("Expected", length(expected) * 3)
)

observed_data <- data.frame(
  Group = rep(c("Metastatic", "Rounded", "Death"), each = length(labels)),
  Category = rep(labels, times = 3),
  Frequency = as.vector(observed),
  Type = rep("Observed", length(as.vector(observed)))
)

combined_data <- rbind(observed_data, expected_data)

# Stacked Bar Plot
ggplot(combined_data, aes(x = Category, y = Frequency, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Group) +
  labs(title = "Observed vs Expected Frequencies",
       x = "Category",
       y = "Frequency") +
  scale_fill_manual(values = c("Observed" = "blue", "Expected" = "red")) +
  theme_minimal()

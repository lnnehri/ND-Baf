library(ggplot2)

#These are finished codes for modeling distributions of the process
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
##THIS SIGMOID FINISHED
####################################################################

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
###########################################################

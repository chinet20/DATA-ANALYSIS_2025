#--------------------start-------------------------------
# Get current working directory
getwd()

#install.packages("PMCMRplus")
#install.packages("coin")
library(PMCMRplus)
library(coin)  # perm.relation
library(dplyr) 

#----------------read dataset--------------------------
data <- read.csv("/Users/ofia/Desktop/Practice 4-20250513/data_for_analysis.csv")
summary(data)

cat("Lipids1", sum(duplicated(data$lipids1)), "\n")
cat("Lipids2 :", sum(duplicated(data$lipids2)), "\n")

# testing for normality of distribution
shapiro.test(data$lipids1)
shapiro.test(data$lipids2)

par(mfrow=c(1,2))
hist(data$lipids1, main = "Lipids1（Histogram of Lipids1）")  
qqnorm(data$lipids1)
qqline(data$lipids1)

# Spearman's correlation test 
spearman_result <- cor.test(data$lipids1, data$lipids2, method="spearman", exact = FALSE)
print(spearman_result)

# data.frame for result
results <- data.frame(
  variable = character(),
  spearman_corr = numeric(),
  s_p_value = numeric(),
  stringsAsFactors = FALSE
)

# variables for analysis
target_vars <- c("lipids2", "lipids3", "lipids4")

# main 
for (var in target_vars) {
  # coinspearman_test
  perm_spearman <- spearman_test(
    data$lipids1 ~ data[[var]],
    distribution = approximate(B = 10000)  # 10000
  )
  
  # add result
  results <- rbind(results, data.frame(
    variable = var,
    spearman_corr = cor(data$lipids1, data[[var]], method = "spearman"),
    s_p_value = pvalue(perm_spearman)
  ))
}

# output result
print(results)

#------visualization of significant results of correlation analysis---------
data <- data[order(data$lipids1),]

plot(data$lipids1, data$lipids2, 
     main = "Lipids1Lipids2（Relationship between Lipids1 and Lipids2）",
     xlab = "Lipids1",
     ylab = "Lipids2",
     pch = 16,
     col = rgb(0,0,1,0.5))
lines(data$lipids1, predict(lm(lipids2 ~ lipids1, data = data)), col = "red", lwd = 2)


lines(data$lipids1, loess(lipids2 ~ lipids1, data = data)$fitted, col = "blue", lwd = 2, lty = 2)


#_____________regression analysis________________ 
df <- data
df <- df[order(df$lipids1),]

#linear regression
model_linear <- lm(lipids1 ~ lipids2, data = df)
summary(model_linear)

#second degree polynomial
model_2 <- lm(lipids1 ~ poly(lipids2, 2, raw = TRUE), data = df)
summary(model_2)

#third degree polynomial
model_3 <- lm(lipids1 ~ poly(lipids2, 3, raw = TRUE), data = df)
summary(model_3)

#exponential dependence 
model_exp <- lm(lipids1 ~ exp(lipids2), data = df)
summary(model_exp)

# log dependence 
model_log <- lm(lipids1 ~ log(lipids2), data = df)
summary(model_log)

#comparison of models
#table of result
rezult <- data.frame(
  model = c("model_linear", "model_2", "model_3", "model_exp", "model_log"),
  BIC_value = c(BIC(model_linear), BIC(model_2), BIC(model_3), BIC(model_exp), BIC(model_log)),
  R_squared = c(summary(model_linear)$r.squared, 
                summary(model_2)$r.squared, 
                summary(model_3)$r.squared, 
                summary(model_exp)$r.squared, 
                summary(model_log)$r.squared)
)

rezult <- rezult[order(rezult$BIC_value),]
print(rezult)

# __________building graphs______________
par(mfrow=c(2,3))

plot(df$lipids2, df$lipids1, 
     main = "（Original Data）", 
     pch = 16, 
     col = rgb(0,0,1,0.5))

plot(df$lipids2, df$lipids1, 
     main = "（Linear Model）", 
     pch = 16, 
     col = rgb(0,0,1,0.5))
lines(df$lipids2, predict(model_linear), col = "red", lwd = 2)

plot(df$lipids2, df$lipids1, 
     main = "（Quadratic Polynomial）", 
     pch = 16, 
     col = rgb(0,0,1,0.5))
lines(df$lipids2, predict(model_2), col = "blue", lwd = 2)


plot(df$lipids2, df$lipids1, 
     main = "（Cubic Polynomial）", 
     pch = 16, 
     col = rgb(0,0,1,0.5))
lines(df$lipids2, predict(model_3), col = "green", lwd = 2)

plot(df$lipids2, df$lipids1, 
     main = "（Exponential Model）", 
     pch = 16, 
     col = rgb(0,0,1,0.5))
lines(df$lipids2, predict(model_exp), col = "purple", lwd = 2)

plot(df$lipids2, df$lipids1, 
     main = "（Logarithmic Model）", 
     pch = 16, 
     col = rgb(0,0,1,0.5))
lines(df$lipids2, predict(model_log), col = "orange", lwd = 2)

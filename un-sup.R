# packages

library(dplyr)
library(corrplot)
library(car)

# Reading the data, basic transformation, EDA

diabetes_data <- read.csv("diabetes.csv", stringsAsFactors = FALSE)
summary(diabetes_data)
head(diabetes_data)
str(diabetes_data)

sapply(diabetes_data, function(x) sum(is.na(x)))
sapply(diabetes_data, class)

par(mfrow = c(3, 9))

target <- diabetes_data
lapply(names(target), function(var) {
  qqnorm(target[[var]], main = '', xlab = '')
  qqline(target[[var]], col = 'lightblue')
  hist(target[[var]], main = var, xlab = '', col = 'lightblue', border = 'black')
  boxplot(target[[var]], ylab = 'Value', col = 'lightblue')
  invisible(NULL)
})


diab_continuous_data <- select(diabetes_data, -c(Outcome)) # continuous variable(s). Excluded initially.
head(diab_continuous_data)

# Perform Shapiro-Wilk test for each continuous variable

shapiro_results <- lapply(names(diab_continuous_data), function(var) {
  test_result <- shapiro.test(diab_continuous_data[[var]])
  c(Statistic = test_result$statistic, P_Value = test_result$p.value)
})

set.seed(123)
shapiro_df <- as.data.frame(do.call(rbind, shapiro_results))
colnames(shapiro_df) <- c("Statistic (W)", "P_Value")
rownames(shapiro_df) <- names(diab_continuous_data)
shapiro_df$TestResult <- ifelse(shapiro_df$P_Value <= 0.05, "Reject", "Fail to Reject")

shapiro_df

# Perform Shapiro-Wilk test for each continuous variable - Log transformed

diab_continuous_data_log <- lapply(diab_continuous_data, function(x) log(x + 1))

shapiro_results_log <- lapply(names(diab_continuous_data_log), function(var) {
  test_result <- shapiro.test(diab_continuous_data_log[[var]])
  c(Statistic = test_result$statistic, P_Value = test_result$p.value)
})

shapiro_df_log <- as.data.frame(do.call(rbind, shapiro_results_log))
colnames(shapiro_df_log) <- c("Statistic (W)", "P_Value")
rownames(shapiro_df_log) <- names(diab_continuous_data_log)

shapiro_df_log$TestResult <- ifelse(shapiro_df_log$P_Value <= 0.05, "Reject", "Fail to Reject")

shapiro_df_log

# Correlation Matrix

cor_matrix <- cor(diab_continuous_data)

par(mar = c(4, 4, 1, 1), cex = 0.9)
corrplot(cor_matrix, method = 'color', 
         col = colorRampPalette(c('blue', 'white', 'red'))(200), 
         addCoef.col = 'black', 
         number.cex = 0.7, 
         format = 'f', 
         digits = 2, 
         tl.col = 'black',  # Color of the variable labels
         tl.cex = 0.7)      # Font size of the labels

# Covariance Matrix

cov_matrix <- cov(diab_continuous_data)
round(cov_matrix, 4)
round(cor_matrix, 4)

# VIF implementation

full_model <- lm(as.formula(paste('Outcome ~', paste(names(diab_continuous_data), collapse = ' + '))),
                 data = diabetes_data)

vif_values <- vif(full_model)

vif_df <- data.frame(
  Variable = names(vif_values),
  vif = vif_values
)

vif_df$vif <- ifelse(vif_df$vif > 10, 'High VIF (Issue)', 'OK')

vif_dif

# PCA

diab_cont_matrix <- as.matrix(diab_continuous_data)
pca <- prcomp(diab_cont_matrix, scale = TRUE)
plot(pca$x[,1], pca$x[,2], xlab = 'PC 1', ylab = 'PC 2')
 
pca.var

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var) * 100, 1)

barplot(pca.var.per, main = 'Scree Plot', xlab = 'Principal Component', ylab = 'Percent Variation',
        names.arg = paste('PC', 1:length(pca.var.per)))

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var / sum(pca.var) * 100, 1)

cumvar <- cumsum(pca.var.per)

residvar <- 100 - cumvar

plot(1:length(residvar), residvar, type = 'b',
     xlab = 'Number of Principal Components',
     ylab = 'Residual Variance (%)',
     main = 'Elbow Plot for Residual Variance')

# pca.data <- data.frame(Sample=rownames(pca$x),
#                        X=pca$x[,1],
#                        Y=pca$x[,2])
# 
# ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
#   geom_text() +
#   xlab(paste('PC1 - ', pca.var.per[1], '%', sep="")) +
#   ylab(paste('PC2 - ', pca.var.per[2], '%', sep="")) +
#   theme_bw() +
#   ggtitle('PCA Graph')

loadings <- pca$rotation
loadings_df <- as.data.frame(loadings)
loadings_df

pca.data

 

# packages

library(dplyr)
library(corrplot)
library(car)
library(factoextra)
library(gridExtra)

# Reading the data, basic transformation, EDA
par(mar = c(2, 2, 2, 2),par(mfrow = c(3,9)))
diabetes_data <- read.csv("diabetes.csv", stringsAsFactors = FALSE)
summary(diabetes_data)
head(diabetes_data)
str(diabetes_data)

sapply(diabetes_data, function(x) sum(is.na(x)))
sapply(diabetes_data, class)



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

par(mar = c(1, 1, 1, 1), mfrow = c(1,1), cex = 0.9)
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

vif_df

# PCA

diab_cont_matrix <- as.matrix(diab_continuous_data) # Convert to matrix
pca <- prcomp(diab_cont_matrix, scale = TRUE) # scale = TRUE for scaling the data
# the output of prcomp is a prcomp object; i.e. a list with several components:
# $dev is the standard deviation of each principal component, $rotation is the matrix of eigenvectors, aka loadings, 
#$center is the center used for scaling, $scale is the scaling used, 
#$x is the matrix of scores, aka the principal components
plot(pca$x[,1], pca$x[,2], xlab = 'PC 1', ylab = 'PC 2') # Plot the first two principal components
print(pca$x[,*])
par(mfrow = c(3,1))
fviz1 <- fviz_pca_var(
  pca,
  col.var = "contrib", # Color by variables
  col.ind = "contrib", # Color by individuals
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Color gradient
  repel = TRUE     # Avoid text overlapping
  )
fviz2 <-fviz_pca_var( # same but with the third and  fourth PC as x and y axes.
  axes = c(3,4), 
  pca,
  col.var = "contrib", # Color by variables
  col.ind = "contrib", # Color by individuals
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Color gradient
  repel = TRUE     # Avoid text overlapping
)
fviz3 <-fviz_pca_var( # same but with the fifth and sixth PC as x and y axes.
  axes = c(5,6), 
  pca,
  col.var = "contrib", # Color by variables
  col.ind = "contrib", # Color by individuals
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Color gradient
  repel = TRUE     # Avoid text overlapping
)
dev.new()
grid.arrange(fviz1, fviz2, fviz3, ncol= 3) # Arrange the plots in a grid
pca.var # Variance of each principal component

pca.var <- pca$sdev^2 # Variance of each principal component
pca.var.per <- round(pca.var/sum(pca.var) * 100, 1) # Percentage of variance of each principal component

barplot(pca.var.per, main = 'Scree Plot', xlab = 'Principal Component', ylab = 'Percent Variation', 
        names.arg = paste('PC', 1:length(pca.var.per))) # Scree plot

pca.var <- pca$sdev^2 # Variance of each principal component
pca.var.per <- round(pca.var / sum(pca.var) * 100, 1) # Percentage of variance of each principal component

cumvar <- cumsum(pca.var.per) # Cumulative variance

residvar <- 100 - cumvar # Residual variance

plot(1:length(residvar), residvar, type = 'b', # Elbow plot
     xlab = 'Number of Principal Components',
     ylab = 'Residual Variance (%)',
     main = 'Elbow Plot for Residual Variance')

# pca.data <- data.frame(Sample=rownames(pca$x),
#                        X=pca$x[,1],
#             http://127.0.0.1:38225/graphics/c0ec00c7-95bb-40c9-aac6-86656b290994.png           Y=pca$x[,2])
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


 


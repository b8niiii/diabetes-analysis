options(warn = -1)

# packages
library(dplyr)
library(corrplot)
library(car)
library(factoextra)
library(ggplot2)
library(gridExtra)
library(plotly)
library(cluster)
library(FactoMineR)

# Reading the data, basic transformation, EDA
diabetes_data <- read.csv('diabetes.csv', stringsAsFactors = FALSE)
summary(diabetes_data)
head(diabetes_data)
str(diabetes_data)

sapply(diabetes_data, function(x) sum(is.na(x)))
sapply(diabetes_data, class)

par(mfrow=c(3,9))
target <- diabetes_data
lapply(names(target), function(var) {
  qqnorm(target[[var]], main = '', xlab = '')
  qqline(target[[var]], col = 'lightblue')
  hist(target[[var]], main = var, xlab = '', col = 'lightblue', border = 'black')
  boxplot(target[[var]], ylab = 'Value', col = 'lightblue')
  invisible(NULL)
})

diab_continuous_data <- select(diabetes_data, -c(Outcome)) # continuous variable(s). Excluded initially.
diab_continuous_data <- scale(diab_continuous_data)
diab_continuous_data <- as.data.frame(diab_continuous_data[apply(diab_continuous_data, 1, function(x) all(abs(x) < 3)), ])

# Perform Shapiro-Wilk test for each continuous variable
set.seed(123)
shapiro_results <- lapply(names(diab_continuous_data), function(var) {
  test_result <- shapiro.test(diab_continuous_data[[var]])
  c(Statistic = test_result$statistic, P_Value = test_result$p.value)
})

shapiro_df <- as.data.frame(do.call(rbind, shapiro_results))
colnames(shapiro_df) <- c("Statistic (W)", "P_Value")
rownames(shapiro_df) <- names(diab_continuous_data)
shapiro_df$TestResult <- ifelse(shapiro_df$P_Value <= 0.05, "Reject", "Fail to Reject")

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
         tl.col = 'black',
         tl.cex = 0.7)

# Covariance Matrix
cov_matrix <- cov(diab_continuous_data)
round(cov_matrix, 4)
round(cor_matrix, 4)

# VIF implementation
full_model <- lm(as.formula(paste('Outcome ~', paste(names(diab_continuous_data), collapse = ' + '))),
                 data = diabetes_data)

vif_values <- vif(full_model)

vif_df <- data.frame(
  variable = names(vif_values),
  vif = vif_values,
  row.names = NULL
)

vif_df$vif_status <- ifelse(vif_df$vif > 10, 'High VIF (Issue)', 'OK')

# PCA
diab_cont_matrix <- as.matrix(diab_continuous_data) # Convert to matrix
pca <- prcomp(diab_cont_matrix, scale = TRUE) # scale = TRUE for scaling the data

# the output of prcomp is a prcomp object; i.e. a list with several components:
# $dev is the standard deviation of each principal component, $rotation is the matrix of eigenvectors, aka loadings, 
# $center is the center used for scaling, $scale is the scaling used, 
# $x is the matrix of scores, aka the principal components

pca.var <- pca$sdev^2 # Variance of each principal component
pca.var.per <- round(pca.var/sum(pca.var) * 100, 1) # Percentage of variance of each principal component

par(mfrow = c(1, 2))
scree_plot <- barplot(pca.var.per, 
                      main = 'Scree Plot', 
                      xlab = 'Principal Component', 
                      ylab = 'Percent Variation', 
                      names.arg = paste('PC', 1:length(pca.var.per)),
                      col = 'lightblue'
)  # Rotate axis labels if needed for better readability

cumvar <- cumsum(pca.var.per) # Cumulative variance per each component
residvar <- 100 - cumvar # Residual variance

elbow_plot <- plot(1:length(residvar), residvar, type = 'n', # Elbow plot
                   xlab = 'Number of Principal Components',
                   ylab = 'Residual Variance (%)',
                   main = 'Elbow Plot - Remaining Residual Variance') # Percentage of unexplained var left after taking #th PC

lines(1:length(residvar), residvar, col = 'black')  
points(1:length(residvar), residvar, col = 'lightblue', pch = 19) 

optimal_PCs_num <- 5
points(optimal_PCs_num, residvar[optimal_PCs_num], col = '#D22B2B', pch = 19, cex = 2)  
text(optimal_PCs_num, residvar[optimal_PCs_num], labels = paste(optimal_PCs_num, 'PCs ~ 80% Variance Explained'), pos = 3, col = '#D22B2B')

par(mfrow = c(1, 1))

rownames(pca$x) <- 1:nrow(pca$x) # Assign row names to the principal components

pca.data <- as.data.frame(pca$x)
pca.data.summary <- data.frame(sample=rownames(pca$x),
                               X=pca$x[,1],
                               Y=pca$x[,2],
                               Z=pca$x[,3],
                               row.names=NULL)

pca.loadings <- data.frame(variable=names(diab_continuous_data),
                           X=pca$rotation[,1],
                           Y=pca$rotation[,2],
                           Z=pca$rotation[,3],
                           row.names=NULL)

# Normalize sample scores to fit within [-1,1] for better visualization
pca.data.summary$X <- pca.data.summary$X / max(abs(pca.data.summary$X))
pca.data.summary$Y <- pca.data.summary$Y / max(abs(pca.data.summary$Y))
pca.data.summary$Z <- pca.data.summary$Z / max(abs(pca.data.summary$Z))

circle <- data.frame(
  x = cos(seq(0, 2 * pi, length.out = 1000)),
  y = sin(seq(0, 2 * pi, length.out = 1000))
)

# PC1 VS PC2
biplot_pc1v2 <- ggplot() +
  # Plot samples as points (scaled between -1 and 1)
  geom_point(data = pca.data.summary, aes(x = X, y = Y), color = 'lightblue', alpha = 0.6) +  
  # Plot variable loadings as arrows
  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = X, yend = Y), 
               arrow = arrow(length = unit(0.2, 'cm')), color = '#D22B2B') +
  geom_text(data = pca.loadings, aes(x = X, y = Y, label = variable), 
            size = 3.5, color = '#D22B2B', vjust = 1.5) +  
  
  # Draw unit circle (optional for reference)
  geom_path(data = circle, aes(x = x, y = y), color = 'lightblue', linetype = 'dashed') +
  # Set axis limits to -1,1 for both samples and variables
  xlab(paste('PC1 - ', pca.var.per[1], '%', sep = '')) +
  ylab(paste('PC2 - ', pca.var.per[2], '%', sep = '')) +
  xlim(-1, 1) + ylim(-1, 1) +  
  coord_fixed() +  # Keep aspect ratio equal
  theme_bw() +
  ggtitle('PCA Biplot (PC1 vs PC2)')

# PC2 VS PC3
biplot_pc2v3 <- ggplot() +
  geom_point(data = pca.data.summary, aes(x = Y, y = Z), color = 'lightblue', alpha = 0.6) +  
  
  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = Y, yend = Z), 
               arrow = arrow(length = unit(0.2, 'cm')), color = '#D22B2B') +
  geom_text(data = pca.loadings, aes(x = Y, y = Z, label = variable), 
            size = 3.5, color = '#D22B2B', vjust = 1.5) +  
  
  geom_path(data = circle, aes(x = x, y = y), color = 'lightblue', linetype = 'dashed') +
  
  xlab(paste('PC2 - ', pca.var.per[2], '%', sep = '')) +
  ylab(paste('PC3 - ', pca.var.per[3], '%', sep = '')) +
  xlim(-1, 1) + ylim(-1, 1) +  
  coord_fixed() + 
  theme_bw() +
  ggtitle('PCA Biplot (PC2 vs PC3)')

# PC1 VS PC3
biplot_pc1v3 <- ggplot() +
  geom_point(data = pca.data.summary, aes(x = X, y = Z), color = 'lightblue', alpha = 0.6) +  
  
  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = X, yend = Z), 
               arrow = arrow(length = unit(0.2, 'cm')), color = '#D22B2B') +
  geom_text(data = pca.loadings, aes(x = X, y = Z, label = variable), 
            size = 3.5, color = '#D22B2B', vjust = 1.5) +  
  
  geom_path(data = circle, aes(x = x, y = y), color = 'lightblue', linetype = 'dashed') +
  
  xlab(paste('PC1 - ', pca.var.per[1], '%', sep = '')) +
  ylab(paste('PC3 - ', pca.var.per[3], '%', sep = '')) +
  xlim(-1, 1) + ylim(-1, 1) +  
  coord_fixed() + 
  theme_bw() +
  ggtitle('PCA Biplot (PC1 vs PC3)')

grid.arrange(biplot_pc1v2, biplot_pc2v3, biplot_pc1v3, ncol = 3)

# PCA 3D plot (PC 1,2,3)
pca_scaled <- pca$x / max(abs(pca$x))
pca.loadings.scaled <- pca$rotation / max(abs(pca$rotation))

pca_plot_threeD <- plot_ly() %>%
  # Add samples (blue points)
  add_trace(x = pca_scaled[,1], 
            y = pca_scaled[,2], 
            z = pca_scaled[,3], 
            type = 'scatter3d', mode = 'markers',
            marker = list(color = 'lightblue', size = 3),
            name = 'Samples') %>%
  
  # Add PC arrows (red vectors)
  add_trace(x = rep(0, nrow(pca.loadings.scaled)), 
            y = rep(0, nrow(pca.loadings.scaled)), 
            z = rep(0, nrow(pca.loadings.scaled)),
            xend = pca.loadings.scaled[,1], 
            yend = pca.loadings.scaled[,2], 
            zend = pca.loadings.scaled[,3], 
            type = 'scatter3d', mode = 'lines',
            line = list(color = '#D22B2B', width = 10),
            name = 'PC Loadings') %>%
  
  # Add variable names at arrow tips
  add_trace(x = pca.loadings.scaled[,1], 
            y = pca.loadings.scaled[,2], 
            z = pca.loadings.scaled[,3], 
            type = 'scatter3d', mode = 'text',
            text = rownames(pca$rotation),
            textposition = 'top center',
            textfont = list(color = '#D22B2B', size = 12),
            name = 'Variables') %>%
  
  # Set layout with proper axis labels
  layout(title = '3D PCA Graph',
         scene = list(xaxis = list(title = paste('PC1 - ', pca.var.per[1], '%')),
                      yaxis = list(title = paste('PC2 - ', pca.var.per[2], '%')),
                      zaxis = list(title = paste('PC3 - ', pca.var.per[3], '%'))))

pca_plot_threeD

# Determine the range of clusters to test (let's try from 1 to 10 clusters)
wcss <- numeric(10)  # To store within-cluster sum of squares
cluster_num <- 1:10
pca_reduced <- pca$x[, 1:optimal_PCs_num] # 5 PCs are sufficient for the var explaination

# Loop to compute K-means and WCSS for different number of clusters
for (i in cluster_num) {
  kmeans_result <- kmeans(pca_reduced, centers = i, nstart = 50)
  wcss[i] <- kmeans_result$tot.withinss  # Store the WCSS value for each number of clusters
}

# Plot Elbow Plot
plot(cluster_num, wcss, type = 'n', 
     xlab = 'Number of Clusters', 
     ylab = 'Within-cluster Sum of Squares (WCSS)', main = 'Elbow Plot - K-means')

lines(cluster_num, wcss, col = 'black')  
points(cluster_num, wcss, col = 'lightblue', pch = 19) 

# Highlight the optimal number of clusters on the plot
optimal_k <- 5
points(optimal_k, wcss[optimal_k], col = '#D22B2B', pch = 19, cex = 2)  # Highlight the optimal k
text(optimal_k, wcss[optimal_k], labels = paste('Optimal k =', optimal_k), pos = 3, col = '#D22B2B')

kmeans_result <- kmeans(pca_reduced, centers = 5, nstart = 50)# Add cluster assignments to pca.data.summary
pca.data.summary$Cluster <- as.factor(kmeans_result$cluster)

cluster_colors <- c('yellow', 'blue', 'green', 'purple', 'orange')

# Update the biplots with proper color mapping
biplot_pc1v2 <- biplot_pc1v2 +
  geom_point(data = pca.data.summary, aes(x = X, y = Y, color = Cluster), alpha = 0.6) +
  scale_color_manual(values = cluster_colors)

biplot_pc2v3 <- biplot_pc2v3 +
  geom_point(data = pca.data.summary, aes(x = Y, y = Z, color = Cluster), alpha = 0.6) +
  scale_color_manual(values = cluster_colors)
  
biplot_pc1v3 <- biplot_pc1v3 +
  geom_point(data = pca.data.summary, aes(x = X, y = Z, color = Cluster), alpha = 0.6) +
  scale_color_manual(values = cluster_colors)
grid.arrange(biplot_pc1v2, biplot_pc2v3, biplot_pc1v3, ncol = 3)

pca_scaled <- as.data.frame(pca_scaled)
pca_scaled$Cluster <- as.factor(kmeans_result$cluster)

# Update the 3D plot with proper color mapping
pca_plot_threeD <- plot_ly() %>%
  add_trace(x = pca_scaled[,1], 
            y = pca_scaled[,2], 
            z = pca_scaled[,3], 
            type = 'scatter3d', mode = 'markers',
            marker = list(color = cluster_colors[as.numeric(pca_scaled$Cluster)], size = 3),  # Convert factor to numeric index
            name = 'Samples') %>%
  
  add_trace(x = rep(0, nrow(pca.loadings.scaled)), 
            y = rep(0, nrow(pca.loadings.scaled)), 
            z = rep(0, nrow(pca.loadings.scaled)),
            xend = pca.loadings.scaled[,1], 
            yend = pca.loadings.scaled[,2], 
            zend = pca.loadings.scaled[,3], 
            type = 'scatter3d', mode = 'lines',
            line = list(color = '#D22B2B', width = 10),
            name = 'PC Loadings') %>%
  
  add_trace(x = pca.loadings.scaled[,1], 
            y = pca.loadings.scaled[,2], 
            z = pca.loadings.scaled[,3], 
            type = 'scatter3d', mode = 'text',
            text = rownames(pca$rotation),
            textposition = 'top center',
            textfont = list(color = '#D22B2B', size = 12),
            name = 'Variables') %>%
  
  layout(title = '3D PCA Graph',
         scene = list(xaxis = list(title = paste('PC1 - ', pca.var.per[1], '%')),
                      yaxis = list(title = paste('PC2 - ', pca.var.per[2], '%')),
                      zaxis = list(title = paste('PC3 - ', pca.var.per[3], '%'))))

pca_plot_threeD

#--------------------------
# # Hierarchical clustering
# # par(mfrow = c(1,1))
# hcpc <- HCPC(pca.data, graph = FALSE) # Hierarchical clustering on principal components
# 
# # #hcpc output is a list with several components: $call is the call to the function
# # # $desc is the description of the clustering, $data is the data used for clustering,
# # # $clust is the clustering results
# hcpc$call$t$nb.clust # Number of clusters
# plot(hcpc, choice = "tree") # Plot the dendrogram
# plot(hcpc, choice = "map") # Plot the map of the clustering
# plot(hcpc, choice = "3D.map") # Plot the 3D map of the clustering
#--------------------------

sil_scores <- numeric()

for (k in 2:length(cluster_num)) {
  # Run k-means clustering
  kmeans_result <- kmeans(pca_reduced, centers = k, nstart = 50)
  
  # Calculate silhouette scores
  sil <- silhouette(kmeans_result$cluster, dist(pca_reduced))
  
  # Store the average silhouette score for this number of clusters
  sil_scores[k-1] <- mean(sil[, 3])  # sil[, 3] contains the silhouette width for each sample
}

plot(2:10, sil_scores, type = "b", 
     xlab = "Number of Clusters", ylab = "Average Silhouette Score", 
     main = "Silhouette Scores for Different Cluster Numbers")

# Use hierarchical methods to compare between methods, number of clusters and linkage types
# Define linkage methods and distance metrics
linkage_methods <- c('complete', 'average', 'single', 'ward.D', 'ward.D2', 'centroid', 'median')
distance_methods <- c('euclidean', 'manhattan', 'canberra', 'binary', 'minkowski')

# Prepare a dataframe to store the results
sil_results_df <- data.frame(method = character(), 
                         distance_method = character(), 
                         linkage_method = character(), 
                         k = integer(), 
                         s_score = numeric(), 
                         stringsAsFactors = FALSE)

# Loop through each combination of distance and linkage methods, and calculate the silhouette score
for (dist_method in distance_methods) {
  for (linkage_method in linkage_methods) {
    for (num_clusters in 2:10) {  # Number of clusters to consider
      # Perform hierarchical clustering
      hc <- hclust(dist(pca_reduced, method = dist_method), method = linkage_method)
      
      # Cut tree to get clusters
      clusters <- cutree(hc, k = num_clusters)
      
      # Calculate silhouette score
      sil_score <- silhouette(clusters, dist(pca_reduced))
      
      # Store the results
      sil_results_df <- rbind(sil_results_df, data.frame(method = 'Hierarchical',
                                                 distance_method = dist_method, 
                                                 linkage_method = linkage_method,
                                                 k = num_clusters, 
                                                 s_score = mean(sil_score[, 3]))) 
    }
  }
}

# View the results dataframe
sil_results_df$s_score <- round(sil_results_df$s_score, 4)
sil_results_df <- sil_results_df[order(-sil_results_df$s_score), ]
sil_results_df


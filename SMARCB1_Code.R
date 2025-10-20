### Setup
rm(list = ls())

# Libraries
library(gt)
library(cluster)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(GGally)
library(factoextra)
library(plotly)
library(randomForest)
library(RColorBrewer)
library(SHAPforxgboost)  
library(xgboost)
library(caret)
library(viridis)
library(umap)            
library(vegan)           

set.seed(43)

### Load data
df <- read.csv("Phenotype_Table.csv", header = TRUE, stringsAsFactors = FALSE)

### Feature selection
col_numbers <- c(
  41, 42, 43, 44, 45, 46, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
  64, 65, 67, 68, 69, 70, 79, 80, 83, 86, 92, 93, 94, 95
)

convert_to_boolean <- function(x) {
  # Convert empties to NA
  x <- ifelse(x == "" | x == "U" | x == " ", NA, x)
  
  # Impute NAs 
  mode_value <- names(sort(table(x), decreasing = TRUE)[1])
  na_count <- sum(is.na(x))
  x[is.na(x)] <- mode_value
  
results <- lapply(df[, col_numbers], convert_to_boolean)

# Imputation summary and boolean matrix
boolean_data <- as.data.frame(lapply(results, function(x) x$values))
imputed_counts <- sapply(results, function(x) x$na_count)
print("Number of imputed values:")
print(imputed_counts)

if (any(is.na(boolean_data))) {
  print("NA values found in the data:")
  print(boolean_data[which(rowSums(is.na(boolean_data)) > 0), ])
  stop("NA values detected.")
}

numeric_data <- as.data.frame(lapply(boolean_data, function(x) ifelse(is.na(x), NA, as.numeric(x))))

if (any(!is.finite(as.matrix(numeric_data)))) {
  print("Non-finite values found in the numeric data. Please check the data conversion steps.")
  print(which(!is.finite(as.matrix(numeric_data)), arr.ind = TRUE))
  stop("Non-finite values detected.")
}

### Labels
row_labels <- sapply(as.numeric(rownames(numeric_data)), function(x) {
  if (x == 1 | x == 2 | x == 3 | x == 5 | x == 6 | x == 7) {
    return("purple")
  } else if (x >= 8 && x <= 12 | x >= 24 && x <= 31) {
    return("green")
  } else if (x >= 13 && x <= 23) {
    return("blue")
  } else {
    return("black")
  }
})
row_labels <- factor(row_labels, levels = c("purple", "green", "blue"))

### Palette and label names
palette <- brewer.pal(n = 3, name = "Set2")
names(palette) <- c("purple", "green", "blue")

label_names <- c(
  purple = "N-Terminal",
  green  = "C-Terminal",
  blue   = "K364del"
)

### PCA
pca_result <- prcomp(numeric_data, scale. = TRUE)
pca_data <- data.frame(pca_result$x, Label = row_labels, Row = rownames(numeric_data))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Label, label = Row)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = palette, labels = c("N-Terminal", "C-Terminal", "K364del")) +
  ggtitle("PCA of SMARCB1 Mutations") +
  theme_minimal() +
  guides(color = guide_legend(title = "Variants", override.aes = list(size = 6, alpha = 1)))
print(pca_plot)

pca_3d_plot <- plot_ly() %>%
  add_trace(
    data = pca_data[pca_data$Label == "purple", ],
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d", mode = "markers+text", text = ~paste(Row),
    textposition = "top right", marker = list(color = palette["purple"]),
    name = "N-Terminal"
  ) %>%
  add_trace(
    data = pca_data[pca_data$Label == "green", ],
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d", mode = "markers+text", text = ~paste(Row),
    textposition = "top right", marker = list(color = palette["green"]),
    name = "C-Terminal"
  ) %>%
  add_trace(
    data = pca_data[pca_data$Label == "blue", ],
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d", mode = "markers+text", text = ~paste(Row),
    textposition = "top right", marker = list(color = palette["blue"]),
    name = "K364del"
  ) %>%
  layout(
    title = "Principal Component Analysis (3D)",
    legend = list(title = list(text = "Variants"), itemsizing = "constant", orientation = "h",
                  x = 0.5, y = -0.1, xanchor = "center"),
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )
print(pca_3d_plot)

# Loadings and variance
loadings <- pca_result$rotation
print("Loadings for each principal component:")
print(loadings)

variance_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
pc1_var <- variance_explained[1]
pc2_var <- variance_explained[2]
pc3_var <- variance_explained[3]

### UMAP (2D and 3D)
umap_result <- umap(numeric_data, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
umap_data <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Label = row_labels,
  Row   = rownames(numeric_data)
)

umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Label, label = Row)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = palette, labels = c("N-Terminal", "C-Terminal", "K364del")) +
  ggtitle("UMAP of SMARCB1 Mutations") +
  theme_minimal() +
  guides(color = guide_legend(title = "Variants", override.aes = list(size = 6, alpha = 1)))
print(umap_plot)

umap_result_3d <- umap(numeric_data, n_neighbors = 15, min_dist = 0.1,
                       metric = "euclidean", n_components = 3)
umap_data_3d <- data.frame(
  UMAP1 = umap_result_3d$layout[, 1],
  UMAP2 = umap_result_3d$layout[, 2],
  UMAP3 = umap_result_3d$layout[, 3],
  Label = row_labels,
  Row   = rownames(numeric_data)
)

umap_3d_plot <- plot_ly() %>%
  add_trace(
    data = umap_data_3d[umap_data_3d$Label == "purple", ],
    x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
    type = "scatter3d", mode = "markers+text", text = ~paste(Row),
    textposition = "top right", marker = list(color = palette["purple"]),
    name = "N-Terminal"
  ) %>%
  add_trace(
    data = umap_data_3d[umap_data_3d$Label == "green", ],
    x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
    type = "scatter3d", mode = "markers+text", text = ~paste(Row),
    textposition = "top right", marker = list(color = palette["green"]),
    name = "C-Terminal"
  ) %>%
  add_trace(
    data = umap_data_3d[umap_data_3d$Label == "blue", ],
    x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
    type = "scatter3d", mode = "markers+text", text = ~paste(Row),
    textposition = "top right", marker = list(color = palette["blue"]),
    name = "K364del"
  ) %>%
  layout(
    title = "3D UMAP of Variants",
    legend = list(title = list(text = "Variants"), itemsizing = "constant", orientation = "h",
                  x = 0.5, y = -0.1, xanchor = "center"),
    scene = list(
      xaxis = list(title = "UMAP1"),
      yaxis = list(title = "UMAP2"),
      zaxis = list(title = "UMAP3")
    )
  )
print(umap_3d_plot)

### Procrustes alignment (UMAP -> PCA) 
pca_data_3d <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  Label = row_labels,
  Row = rownames(numeric_data)
)

common_rows <- intersect(pca_data_3d$Row, umap_data_3d$Row)
pca_data_3d  <- pca_data_3d [match(common_rows, pca_data_3d$Row), ]
umap_data_3d <- umap_data_3d[match(common_rows, umap_data_3d$Row), ]

pca_coords  <- as.matrix(pca_data_3d[, c("PC1", "PC2", "PC3")])
umap_coords <- as.matrix(umap_data_3d[, c("UMAP1", "UMAP2", "UMAP3")])

proc <- procrustes(pca_coords, umap_coords, scale = TRUE)
aligned_umap <- proc$Yrot

umap_data_3d$UMAP1 <- aligned_umap[, 1]
umap_data_3d$UMAP2 <- aligned_umap[, 2]
umap_data_3d$UMAP3 <- aligned_umap[, 3]

lighter_palette <- sapply(palette, function(col) adjustcolor(col, alpha.f = 0.5))

p <- plot_ly() %>%
  add_trace(
    data = umap_data_3d, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
    type = "scatter3d", mode = "markers",
    color = ~Label, colors = lighter_palette,
    marker = list(size = 5, symbol = "square"),
    name = "UMAP"
  ) %>%
  add_trace(
    data = pca_data_3d, x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d", mode = "markers",
    color = ~Label, colors = palette,
    marker = list(size = 5, symbol = "circle"),
    name = "PCA"
  )

for (i in seq_len(nrow(umap_data_3d))) {
  current_label <- umap_data_3d$Label[i]
  row_name <- umap_data_3d$Row[i]
  
  pca_match <- pca_data_3d[pca_data_3d$Row == row_name, ]
  
  x_umap <- umap_data_3d$UMAP1[i]
  y_umap <- umap_data_3d$UMAP2[i]
  z_umap <- umap_data_3d$UMAP3[i]
  
  x_pca <- pca_match$PC1
  y_pca <- pca_match$PC2
  z_pca <- pca_match$PC3
  
  line_color <- palette[as.character(current_label)]
  
  p <- p %>% add_trace(
    x = c(x_umap, x_pca),
    y = c(y_umap, y_pca),
    z = c(z_umap, z_pca),
    type = "scatter3d", mode = "lines",
    line = list(color = line_color, width = 2),
    showlegend = FALSE
  )
}

p <- p %>%
  layout(
    title = "PCA/UMAP Plot",
    legend = list(title = list(text = "Variants"), orientation = "h",
                  x = 0.5, y = -0.1, xanchor = "center"),
    scene = list(
      xaxis = list(title = "Dimension 1"),
      yaxis = list(title = "Dimension 2"),
      zaxis = list(title = "Dimension 3")
    )
  )
print(p)

### XGBoost
label_numeric <- as.numeric(row_labels) - 1
dtrain <- xgb.DMatrix(data = as.matrix(numeric_data), label = label_numeric)

params <- list(
  booster = "gbtree",
  objective = "binary:logistic",
  eta = 0.1,
  max_depth = 6,
  eval_metric = "logloss"
)

xgb_model <- xgb.train(params = params, data = dtrain, nrounds = 100)

# SHAP analysis
shap_values <- shap.values(xgb_model, X_train = as.matrix(numeric_data))
shap_long <- shap.prep(xgb_model, X_train = as.matrix(numeric_data))

shap.plot.summary(shap_long)

shap_importance <- shap_long %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = mean(abs(value)), .groups = "drop") %>%
  arrange(desc(mean_abs_shap))

cat("Feature importance (mean absolute SHAP values):\n")
print(shap_importance)

palette_vars <- viridis(length(unique(shap_long$variable)), option = "D")
top_features <- shap_importance$variable[1:5]

superimposed_beeswarm_plot <- ggplot(
  shap_long %>% filter(variable %in% top_features),
  aes(x = value, y = reorder(variable, value, FUN = median), color = variable)
) +
  geom_jitter(width = 0.2, height = 0.1, size = 1.5, alpha = 0.6) +
  labs(title = "Superimposed SHAP Bee-Swarm Plot for Top Features",
       x = "SHAP Value", y = "Feature") +
  theme_minimal() +
  scale_color_viridis_d(name = "Features", option = "D") +
  theme(legend.position = "right")
print(superimposed_beeswarm_plot)

### Elastic Net

# Keep only purple/green for binary regression
en_rows <- row_labels %in% c("purple", "green")
X_en <- as.matrix(numeric_data[en_rows, , drop = FALSE])
y_en <- as.numeric(row_labels[en_rows] == "purple")  # 1 = purple, 0 = green

cat("Elastic Net: X rows/cols =", nrow(X_en), "/", ncol(X_en), " | y length =", length(y_en), "\n")

# Standardise predictors (glmnet standardises internally by default; we mirror your approach)
X_en <- scale(X_en)

set.seed(43)
alpha_values <- seq(0, 1, by = 0.1)
best_alpha <- NA_real_
best_model <- NULL
best_lambda <- NA_real_
best_deviance <- Inf

# Cross-validate across alpha values
for (alpha in alpha_values) {
  m <- cv.glmnet(X_en, y_en, alpha = alpha, family = "binomial")
  this_min <- min(m$cvm)
  if (this_min < best_deviance) {
    best_deviance <- this_min
    best_alpha <- alpha
    best_model <- m
    best_lambda <- m$lambda.min
  }
}

# Coefficients at optimal lambda/alpha
elastic_net_coefficients <- as.matrix(coef(best_model, s = best_lambda))
cat(paste0("Best alpha: ", best_alpha, "\n"))
cat("Elastic Net Coefficients (including intercept):\n")
print(elastic_net_coefficients)

# Build tidy coefficient table (exclude intercept, keep non-zero)
coef_df <- data.frame(
  Feature = rownames(elastic_net_coefficients),
  Coefficient = as.vector(elastic_net_coefficients),
  row.names = NULL,
  check.names = FALSE
)

coef_df <- subset(coef_df, Feature != "(Intercept)")
selected_features <- subset(coef_df, Coefficient != 0)

# Pretty names
selected_features$Feature <- gsub("\\.", " ", selected_features$Feature)

# Order by absolute magnitude
selected_features <- selected_features[order(abs(selected_features$Coefficient), decreasing = TRUE), ]

# attach imputation count
# selected_features$ImputedCount <- imputed_counts[selected_features$Feature]

# Plot
if (nrow(selected_features) > 0) {
  print("Selected Features by Elastic Net:")
  print(selected_features)
  
  color_palette <- RColorBrewer::brewer.pal(n = 3, name = "Set2")
  
  ggplot(selected_features,
         aes(x = reorder(Feature, abs(Coefficient)),
             y = Coefficient,
             fill = Coefficient > 0)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = c("TRUE" = color_palette[1], "FALSE" = color_palette[2]),
      labels = c("FALSE" = "C-Terminal", "TRUE" = "N-Terminal"),
      name   = "Variant Location"
    ) +
    coord_flip() +
    ggtitle("Feature Importance of C-Terminal and N-Terminal Variants of SMARCB1 (Elastic Net)") +
    xlab("Phenotype") +
    ylab("Coefficient") +
    theme_minimal()
} else {
  message("No features were selected by Elastic Net.")
}

### Random Forest

pca_data <- data.frame(pca_result$x, Label = row_labels, Row = rownames(numeric_data))

set.seed(43)
rf_model <- randomForest(pca_data[, 1:3], row_labels, ntree = 100, mtry = 2, importance = TRUE)
print(rf_model)
plot(rf_model, main = "OOB Error Rate by Number of Trees")
varImpPlot(rf_model, main = "Variable Importance")

# Refit with proximity for MDS
set.seed(42)
rf_model <- randomForest(pca_data[, 1:3], row_labels,
                         ntree = 100, mtry = 2, importance = TRUE, proximity = TRUE)
MDSplot(rf_model, row_labels, palette = c("purple", "green"), pch = 19)
title("MDS Plot of Random Forest Proximity")
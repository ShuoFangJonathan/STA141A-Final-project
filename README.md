# STA141A-Final-project
```{r}
# install packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyverse)
library(forcats)
library(caret)
library(gridExtra)
library(xgboost)
library(pROC)
```


```{r}
session <- list()
base_path <- "/Users/jonathanfang/Desktop/STA 141A Final Report/sessions/"
for(i in 1:18) {
  file_path <- paste0(base_path, "session", i, ".rds")
  session[[i]] <- readRDS(file_path)
  print(session[[i]]$mouse_name)
  print(session[[i]]$date_exp)
}

session_summary <- map_df(1:18, function(i) {
  session_data <- session[[i]]
  tibble(
    SessionID = i,
    MouseName = session_data$mouse_name,
    DateExp = session_data$date_exp,
    NumberOfNeurons = length(unique(session_data$brain_area)),
    NumberOfTrials = length(session_data$spks),
    AverageSuccessRate = mean(session_data$feedback_type == 1, na.rm = TRUE)
  )
})
print(session_summary)
```
```{r}
names(session[[8]])
dim(session[[8]]$spks[[1]]) 
length(session[[1]]$brain_area)
session[[8]]$spks[[1]][6,]
session[[8]]$spks[[1]][6,3] 
session[[8]]$brain_area[6]
```
```{r}
get_trail_data <- function(session_id, trail_id){
  spikes <- session[[session_id]]$spks[[trail_id]]
  if (any(is.na(spikes))){
    disp("value missing")
  }

  #trail_tibble <- as_tibble(spikes) %>% set_names(binename) %>%  add_column("brain_area" = session[[session_id]]$brain_area ) %>% group_by(brain_area) %>% summarize( "sum_spikes" =across(everything(),sum),.groups = "drop") 
  trail_tibble <- tibble("neuron_spike" = rowSums(spikes))  %>%  add_column("brain_area" = session[[session_id]]$brain_area ) %>% group_by(brain_area) %>% summarize( region_sum_spike = sum(neuron_spike), region_count = n(),region_mean_spike = mean(neuron_spike)) 
  trail_tibble  = trail_tibble%>% add_column("trail_id" = trail_id) %>% add_column("contrast_left"= session[[session_id]]$contrast_left[trail_id]) %>% add_column("contrast_right"= session[[session_id]]$contrast_right[trail_id]) %>% add_column("feedback_type"= session[[session_id]]$feedback_type[trail_id])
  trail_tibble
}
trail_tibble_1_2 <- get_trail_data(1,2)
trail_tibble_1_2

get_session_data <- function(session_id){
  n_trail <- length(session[[session_id]]$spks)
  trail_list <- list()
  for (trail_id in 1:n_trail){
    trail_tibble <- get_trail_data(session_id,trail_id)
    trail_list[[trail_id]] <- trail_tibble
  }
  session_tibble <- do.call(rbind, trail_list)
  session_tibble <- session_tibble %>% add_column("mouse_name" = session[[session_id]]$mouse_name) %>% add_column("date_exp" = session[[session_id]]$date_exp) %>% add_column("session_id" = session_id) 
  session_tibble
}

session_1 <- get_session_data(1)
head(session_1)

session_list = list()
for (session_id in 1: 18){
  session_list[[session_id]] <- get_session_data(session_id)
}
full_tibble <- do.call(rbind, session_list)
full_tibble$success <- full_tibble$feedback_type == 1
full_tibble$success <- as.numeric(full_tibble$success)
full_tibble$contrast_diff <- abs(full_tibble$contrast_left-full_tibble$contrast_right)
```
```{r}
full_tibble %>% filter (trail_id==1) %>% group_by(session_id) %>% summarise(sum(region_count))
full_tibble %>% group_by(session_id) %>% summarise(unique_area = n_distinct(brain_area))
average_spike <-full_tibble %>% group_by( session_id, trail_id) %>% mutate(mean_spike = sum(region_sum_spike)/sum(region_count))
average_spike %>% group_by(session_id) %>% summarise(mean_session_spike = mean(mean_spike))
full_tibble$session_id <- factor(full_tibble$session_id)
ggplot(full_tibble, aes(x = session_id, y = brain_area)) +
  geom_point() +
  labs(x = "Session ID", y = "Brain Area") +
  theme_minimal()
full_functional_tibble %>% group_by(session_id) %>% summarize(success_rate = mean(success, na.rm = TRUE))
full_functional_tibble %>% group_by(mouse_name) %>% summarize(success_rate = mean(success, na.rm = TRUE))

full_functional_tibble %>% group_by(contrast_diff) %>% count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

full_functional_tibble %>% group_by(contrast_diff) %>% summarize(success_rate = mean(success, na.rm = TRUE))

counts_df <- full_functional_tibble[c('mouse_name', 'contrast_diff')]
counts_df$contrast_diff <- as.factor(counts_df$contrast_diff)
counts <- table(counts_df)

percentages <- prop.table(counts, margin = 1)
percentages


full_functional_tibble$contrast_diff <- as.factor(full_functional_tibble$contrast_diff)
anova_contrast_mouse <- aov(success ~ mouse_name * contrast_diff, data = full_functional_tibble)
summary(anova_contrast_mouse)


```

# visulize
```{r}
full_functional_tibble$trial_group = cut(
  full_functional_tibble$trial_id,
  breaks = seq(0, max(full_functional_tibble$trial_id, na.rm = TRUE), by = 25),
  include.lowest = TRUE
)


success_rate <- aggregate(success ~ session_id + trial_group, data = full_functional_tibble, FUN = mean)
ggplot(success_rate, aes(x = trial_group, y = success)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~session_id, ncol=3) +
  theme_bw()

success_rate <- aggregate(success ~ mouse_name + trial_group, data = full_functional_tibble, FUN = mean)
ggplot(success_rate, aes(x = trial_group, y = success)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~mouse_name) +
  theme_bw()

```
```{r}
col_names <-names(full_functional_tibble)
region_sum_subset <- col_names[grep("^region_sum", col_names)]
region_mean_subset <- col_names[grep("^region_mean", col_names)]

average_spike <- full_tibble %>% group_by( session_id,trail_id) %>% summarise(mean_spike = sum(region_sum_spike)/sum(region_count))

average_spike$mouse_name <- full_functional_tibble$mouse_name
average_spike$contrast_diff <- full_functional_tibble$contrast_diff
average_spike$success <- full_functional_tibble$success

ggplot(average_spike, aes(x = trail_id, y = mean_spike)) + 
  geom_line()+
  geom_smooth(method = "loess")+

  facet_wrap(~session_id)

ggplot(average_spike, aes(x = trail_id, y = mean_spike)) + 
  geom_line()+
  geom_smooth(method = "loess")+

  facet_wrap(~mouse_name)
```
```{r}
## PCA and visualize the 2D results
features = full_functional_tibble[,1:40]
scaled_features <- scale(features)
pca_result <- prcomp(scaled_features)
pc_df <- as.data.frame(pca_result$x)
pc_df$session_id <- full_functional_tibble$session_id
pc_df$mouse_name <- full_functional_tibble$mouse_name

p1 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = session_id)) +
  geom_point(size = 1, alpha = 0.4, shape = 16) + 
  labs(title = "PCA: PC1 vs PC2.Session ID") +
  theme_classic() +
  theme(legend.position = "right", aspect.ratio = 1)
p2 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = mouse_name)) +
  geom_point(size = 1, alpha = 0.4, shape = 16) + 
  labs(title = "PCA: PC1 vs PC2.Mouse Name") +
  theme_classic() +
  theme(legend.position = "right", aspect.ratio = 1)
grid.arrange(p1, p2, ncol = 2)
```
##  Data integration


```{r}
colnames(full_functional_tibble)
predictive_feature <- c("session_id", "trial_id", "contrast_right", "contrast_left", "contrast_diff","binename")
head(full_functional_tibble[predictive_feature])
```
```{r}
auroc <- roc(test_label, predictions)
auroc
```





```{r}
# split model for session 15
set.seed(123) # for reproducibility
session_15_row <- which(full_functional_tibble$session_id==15)
testIndex <- sample(session_15_row, 50, replace = FALSE)
trainIndex <- 1:nrow(full_functional_tibble)
trainIndex <- trainIndex[!(trainIndex %in% testIndex)]

train_df <- predictive_dat[trainIndex, ]
train_X <- X[trainIndex,]
test_df <- predictive_dat[-trainIndex, ]
test_X <- X[-trainIndex,]

train_label <- label[trainIndex]
test_label <- label[-trainIndex]

xgb_model <- xgboost(data = train_X, label = train_label, objective = "binary:logistic", nrounds=10)
predictions <- predict(xgb_model, newdata = test_X)
predicted_labels <- as.numeric(ifelse(predictions > 0.5, 1, 0))
accuracy <- mean(predicted_labels == test_label)
accuracy
conf_matrix <- confusionMatrix(as.factor(predicted_labels), as.factor(test_label))
conf_matrix$table
auroc <- roc(test_label, predictions)
auroc
```



```{r}
# split model for session 2
set.seed(123) # for reproducibility
session_2_row <- which(full_functional_tibble$session_id==2)
testIndex <- sample(session_2_row, 50, replace = FALSE)
trainIndex <- 1:nrow(full_functional_tibble)
trainIndex <- trainIndex[!(trainIndex %in% testIndex)]

train_df <- predictive_dat[trainIndex, ]
train_X <- X[trainIndex,]
test_df <- predictive_dat[-trainIndex, ]
test_X <- X[-trainIndex,]

train_label <- label[trainIndex]
test_label <- label[-trainIndex]
xgb_model <- xgboost(data = train_X, label = train_label, objective = "binary:logistic", nrounds=10)
predictions <- predict(xgb_model, newdata = test_X)
predicted_labels <- as.numeric(ifelse(predictions > 0.5, 1, 0))
accuracy <- mean(predicted_labels == test_label)
accuracy
conf_matrix <- confusionMatrix(as.factor(predicted_labels), as.factor(test_label))
conf_matrix$table
auroc <- roc(test_label, predictions)
auroc
```



```{r}
results <- data.frame(
  Model = c("predictive_dat", "session_15_row", "session_2_row"),
  Accuracy = c(0.7224, 0.82, 0.62),
  AUC = c(0.7254, 0.813, 0.6929)
)

accuracy_plot <- ggplot(results, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity") +
  ylim(0, 1) +
  labs(title = "Model Accuracy", y = "Accuracy", x = "") +
  theme_minimal()

auc_plot <- ggplot(results, aes(x = Model, y = AUC, fill = Model)) +
  geom_bar(stat = "identity") +
  ylim(0, 1) +
  labs(title = "Model AUC", y = "AUC", x = "") +
  theme_minimal()
conf_matrix1 <- matrix(c(60, 45, 237, 674), nrow = 2)
conf_matrix2 <- matrix(c(2, 2, 7, 39), nrow = 2)
conf_matrix3 <- matrix(c(4, 2, 17, 27), nrow = 2)

conf_matrices <- list(conf_matrix1, conf_matrix2, conf_matrix3)
plot_conf_matrix <- function(cm, model_name) {
  cm <- as.table(cm)
  ggplot(data = as.data.frame(cm), aes(x = Var1, y = Var2)) +
    geom_tile(aes(fill = Freq), colour = "white") +
    geom_text(aes(label = Freq), vjust = 1) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = paste("Confusion Matrix:", model_name),
         x = "Predicted Label", y = "Actual Label") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

conf_matrix_plots <- lapply(seq_along(conf_matrices), function(i) {
  plot_conf_matrix(conf_matrices[[i]], paste("Model", i))
})

library(gridExtra)
grid.arrange(accuracy_plot, auc_plot, conf_matrix_plots[[1]], conf_matrix_plots[[2]], conf_matrix_plots[[3]], ncol = 2)
```

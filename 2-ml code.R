library(tidyverse)
library(dplyr) 
library(MASS)

library(mlbench)

data("BreastCancer")


# Check first few rows
head(BreastCancer)


names(BreastCancer)
head(BreastCancer)
dim(BreastCancer)
summary(BreastCancer)

bc <- BreastCancer[, -1]
bc <- bc %>%
  mutate(across(-Class, ~ as.numeric(as.character(.))))
glimpse(bc)

# Loop over each variable in bc except 'Class'
for (var in names(bc)[names(bc) != "Class"]) {
  # Create a table of counts for the current variable
  counts <- table(bc[[var]])
  
  # Plot bar chart using base R
  barplot(counts,
          main = paste("Bar Chart of", var),
          xlab = var,
          ylab = "Frequency",
          col = "skyblue",
          border = "black")
}

# Count missing values in each column
colSums(is.na(bc))
bc$Bare.nuclei[is.na(bc$Bare.nuclei)]<-mean(bc$Bare.nuclei,na.rm=TRUE)

# Set seed and split
set.seed(123)
n <- nrow(bc)
train_indices <- sample(seq_len(n), size = 0.7 * n)
train_data <- bc[train_indices, ]
test_data  <- bc[-train_indices, ]

# Separate features and target
X_train <- train_data[, setdiff(names(train_data), "Class")]
y_train <- train_data$Class

X_test <- test_data[, setdiff(names(test_data), "Class")]
y_test <- test_data$Class

# Standardize features
X_train_scaled <- as.data.frame(scale(X_train))
X_test_scaled <- as.data.frame(scale(X_test,
                                     center = attr(scale(X_train), "scaled:center"),
                                     scale  = attr(scale(X_train), "scaled:scale")))

# Recombine with Class
train_data_std <- cbind(X_train_scaled, Class = y_train)
test_data_std  <- cbind(X_test_scaled,  Class = y_test)



# LDA
lda_model <- lda(Class ~ ., data = train_data_std)
lda_model
lda_pred <- predict(lda_model, test_data_std)$class

lda_cm <- table(Predicted = lda_pred, Actual = test_data_std$Class)
print(lda_cm)

# QDA
qda_model <- qda(Class ~ ., data = train_data_std)
qda_pred <- predict(qda_model, test_data_std)$class

qda_cm <- table(Predicted = qda_pred, Actual = test_data_std$Class)
print(qda_cm)

# Logistic Regression
log_model <- glm(Class ~ ., data = train_data_std, family = binomial)
log_probs <- predict(log_model, newdata = test_data_std, type = "response")
# Convert probs to class labels (adjust threshold as needed)
log_pred <- ifelse(log_probs >= 0.5, levels(y_train)[2], levels(y_train)[1])
log_pred <- factor(log_pred, levels = levels(y_train))
log_cm <- table(Predicted = log_pred, Actual = test_data_std$Class)
print(log_cm)


#defining a function for model evaluation
evaluate_model <- function(predicted, actual) {
  cm <- table(Predicted = predicted, Actual = actual)
  
  accuracy <- sum(diag(cm)) / sum(cm)
  precision <- diag(cm) / rowSums(cm)
  recall <- diag(cm) / colSums(cm)
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  
  # Macro-averaged metrics: average across classes
  macro_precision <- mean(precision)
  macro_recall <- mean(recall)
  macro_f1 <- mean(f1)
  
  # Return as named vector or data frame
  return(c(
    Accuracy = round(accuracy, 4),
    Precision = round(macro_precision, 4),
    Recall = round(macro_recall, 4),
    F1_Score = round(macro_f1, 4)
  ))
}
lda_results <- evaluate_model(lda_pred, test_data_std$Class)
lda_results
qda_results <- evaluate_model(qda_pred, test_data_std$Class)
qda_results 
log_results <- evaluate_model(log_pred, test_data_std$Class)
log_results
comparison_table <- rbind(
  LDA = lda_results,
  QDA = qda_results,
  Logistic = log_results
)

print(comparison_table)
library(pROC)
roc_lda <- roc(test_data_std$Class, lda_probs)
roc_qda <- roc(test_data_std$Class, qda_probs)
roc_log <- roc(test_data_std$Class, log_probs)

# --- Plot ---
plot(roc_log, col="blue", main="ROC Curves for LDA, QDA, Logistic")
lines(roc_lda, col="green")
lines(roc_qda, col="red")
legend("bottomright", legend=c("Logistic", "LDA", "QDA"),
       col=c("blue", "green", "red"), lwd=2)



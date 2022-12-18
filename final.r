# load packages 
library(pls)
library(tidyverse)
library(caret)

# Load the data
  data("Boston", package = "MASS") # or just simply call MASS package from the beginning... 

# Split the data into training and test set
set.seed(123)
training.samples <- Boston$medv %>%
  createDataPartition(p = 0.8, list = FALSE)

train.data  <- Boston[training.samples, ]
test.data <- Boston[-training.samples, ]

# 1.Principal Component Regression (PCR) ---------------------

# Build the model on training set (method = pcr) 
set.seed(123)
model.pcr <- train(
  medv~., data = train.data, method = "pcr",
  scale = TRUE,
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)

# Plot model RMSE vs different values of components
plot(model.pcr)

# Print the best tuning parameter ncomp that minimize the cross-validation error(->RMSE) 
model.pcr$bestTune

# Make predictions
pred.pcr <- model.pcr %>% predict(test.data)
# Model performance metrics
data.frame(
  RMSE = caret::RMSE(pred.pcr, test.data$medv),
  Rsquare = caret::R2(pred.pcr, test.data$medv)
)


# 2. Partial Least Squares Regression (PLSR) ------------------------------------------------

# Build the model on training set (method = pls) 
set.seed(123)
model.pls <- train(
  medv~., data = train.data, method = "pls",
  scale = TRUE,
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)
# Plot model RMSE vs different values of components
plot(model.pls)
# Print the best tuning parameter ncomp that
# minimize the cross-validation error, RMSE
model.pls$bestTune

# Summarize the final model
summary(model.pls$finalModel)

# Make predictions
pred.pls <- model.pls %>% predict(test.data)
# Model performance metrics
data.frame(
  RMSE = caret::RMSE(pred.pls, test.data$medv),
  Rsquare = caret::R2(pred.pls, test.data$medv)
)


# 3. LASSO Regression (LR) ------------------------------------------------

# Build the model on training set (method = lasso) 
set.seed(123)
model.lasso <- train(
  medv~., data = train.data, method = "glmnet",
  scale = TRUE,
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)
# Plot model RMSE vs different values of components
plot(model.lasso)
# Print the best tuning parameter ncomp that
# minimize the cross-validation error, RMSE
model.lasso$bestTune

# Summarize the final model
summary(model.lasso$finalModel)

# Make predictions
pred.lasso <- model.lasso %>% predict(test.data)
# Model performance metrics
data.frame(
  RMSE = caret::RMSE(pred.lasso, test.data$medv),
  Rsquare = caret::R2(pred.lasso, test.data$medv)
)


# 4. Random Forest (RF) -----------------------------------

# Build the model on training set (method = pls) 
set.seed(123)
model.rf <- train(
  medv~., data = train.data, method = "rf",
  scale = TRUE,
  trControl = trainControl("cv", number = 100),
  tuneLength = 10
)
# Plot model RMSE vs different values of components
plot(model.rf)
# Print the best tuning parameter ncomp that
# minimize the cross-validation error, RMSE
model.rf$bestTune

# Summarize the final model
summary(model.rf$finalModel)

# Make predictions
pred.rf <- model.rf %>% predict(test.data)
# Model performance metrics
data.frame(
  RMSE = caret::RMSE(pred.pls, test.data$medv),
  Rsquare = caret::R2(pred.pls, test.data$medv)
)

plot(model.pcr)
par(new=TRUE)
plot(model.pls)
lines(model.rf)

## to sum up: PLS the best (supervised method), RF the worst result (time consuming, low accuracy)

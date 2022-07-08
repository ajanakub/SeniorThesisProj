## LassoPractice.R - Practicing Lasso.
## Bukola Ajanaku
## February 3, 2021
## Practicing lasso using blog tutorial. Following tutorial, combination of
## both tutorial from blog and tutorial from youtube. This following tutorial
## is freestyled by me.

library(glmnet)

# creating the variables from the mtcars innate R dataset
responseVar <- mtcars$hp
predictorVar <- data.matrix(mtcars[, c("mpg", "wt", "drat", "qsec")])

# Creates optimal lambda values using a "k-fold" cross-validation:
# k-fold validation breaks the data set into a k amount of groups
# (aka folds). Then it tests each fold against the model data set
# to develop the mean squared error (MSE) to determine how well the given
# model would work on a dataset that it hasn't previously seen. Thus,
# this step works towards recreating a model by identifying lambda
# values that allow for the lowest MSE. The cv.glmnet uses k = 10 folds.
# Alpha is from 0 to 1 which determines how tightly the variables should
# be correlated. Elastic net (a = 1) is looser and reduces the impact of different
# features without eliminating them whereas Ridge (a = 0) tries to minimize the
# impact of features that are deemed not important.
crossvalid <- cv.glmnet(predictorVar, responseVar, alpha = 1)

# Now, we want to grab only the lambda values that minimize the test MSE
# values for the model we are creating.
mainlambda <- crossvalid$lambda.min
mainlambda

# Now, it's time to create the model. Any variabl that does not have a coeff
# was shrunken all the way to zero because it didn't hold enough importance.
model <- glmnet(predictorVar, responseVar, alpha = 1, lambda = mainlambda)
coef(model)

# Now, we can use the lassoed model to make predictions.
responseVar_predicted <- predict(model, s = mainlambda, newx = predictorVar)

# Finding R-squared. First, we'd need to calculate the sum of squares total,
# sum of squares error, to get the r-squared.
SST <- sum((responseVar - mean(responseVar))^2)
SSE <- sum((responseVar_predicted-responseVar)^2)

RSQ <- 1 - (SSE/SST)
RSQ

# The RSQ value shows the percentage of the data that can be explained by the
# model we came up with!

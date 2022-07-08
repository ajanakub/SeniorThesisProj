## LassOP.R - Practicing Lasso from https://www.youtube.com/watch?v=ctmNq7FgbvI
## Bukola Ajanaku
## February 3, 2021

library(glmnet)

responseVar <- mtcars$hp
predictorVar <- data.matrix(mtcars[, c("mpg", "wt", "drat", "qsec")])

# Ridge Regression example
alpha0.fit <- cv.glmnet(predictorVar, responseVar, type.measure = "mse",
alpha = 0, family = "gaussian")

alpha0.predicted <- predict(alpha0.fit, s = alpha0.fit$lambda.min, newx = predictorVar)
mean((responseVar - alpha0.predicted)^2)

# Lasso Regression example
alpha1.fit <- cv.glmnet(predictorVar, responseVar, type.measure = "mse",
alpha = 1, family = "gaussian")

alpha1.predicted <- predict(alpha1.fit, s = alpha1.fit$lambda.min, newx = predictorVar)
mean((responseVar - alpha1.predicted)^2)


# similar to spar
test_that("Results has right class", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  y <- rnorm(10)
  spar_res <- spar.cv(x,y, rp=rp_gaussian(),model = spar_glm())
  expect_equal(class(spar_res),"spar.cv")
})

test_that("Coef returns vector of correct length", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar.cv(x,y, model = spar_glm())
  sparcoef <- coef(spar_res)
  expect_equal(length(sparcoef$beta),30)
})

test_that("Coef is more sparse for 1se rule", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- rnorm(10)
  spar_res <- spar.cv(x,y,screencoef = screen_glmnet(),
                      model = spar_glm())
  sparcoef <- coef(spar_res,opt_par = "best")
  sparcoef2 <- coef(spar_res,opt_par = "1se")
  expect_equal(all(which(sparcoef$beta==0) %in% which(sparcoef2$beta==0)),TRUE)
})

test_that("Validated nu values are same as the ones for initial SPAR fit", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  y <- rnorm(10)
  spar_res <- spar.cv(x,y,screencoef = screen_glmnet(),
                      nummods=c(10,15), model = spar_glm())
  expect_equal(unique(spar_res$val_sum$nu),as.numeric(spar_res$nus))
})

test_that("Validated nummod values are same as the ones for initial SPAR fit", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  y <- rnorm(10)
  spar_res <- spar.cv(x,y,screencoef = screen_glmnet(),
                      nummods=c(10,15), model = spar_glm())
  expect_equal(unique(spar_res$val_sum$nummod),as.numeric(spar_res$nummods))
})

test_that("Columns with zero sd get coefficient 0", {
  x <- example_data$x
  x[,c(1,11,111)] <- 2
  y <- example_data$y
  spar_res <- spar.cv(x, y, screencoef = screen_glmnet(),
                      measure = "mae", model = spar_glm())
  sparcoef <- coef(spar_res)
  expect_equal(sparcoef$beta[c(1,11,111)],c(0,0,0))
})

# Tests expecting errors

test_that("Get errors for input x not data.frame or matrix", {
  x <- list("1"=1:10,"2"=(-1)^(1:12),"3"=rnorm(12),
            "4"=rnorm(12),"5"=runif(12),"6"=runif(12))
  y <- rnorm(6)
  expect_error(spar.cv(x,y,nfolds=2, model = spar_glm()))
})

test_that("Get errors for input x non-numeric data.frame", {
  x <- data.frame(matrix(rnorm(300), ncol = 30))
  x$X1[1] <- "a"
  y <- rnorm(10)
  expect_error(spar.cv(x,y))
})

test_that("Get errors for categorical input y", {
  x <- matrix(rnorm(300), ncol = 30)
  y <- factor(rep(c("a","b"),5))
  expect_error(spar.cv(x,y))
})

test_that("Get errors for prediction when xnew has wrong dimensions", {
  x <- example_data$x
  y <- example_data$y
  spar_res <- spar.cv(x,y,model = spar_glm())
  xnew <- example_data$xtest
  expect_error(predict(spar_res,xnew=xnew[,-1]))
})

test_that("Get errors for classification validation measure for non-binomial family", {
  x <- example_data$x
  y <- example_data$y
  expect_error(spar.cv(x,y,measure = "1-auc",model = spar_glm()))
})




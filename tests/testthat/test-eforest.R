

# Load data ---------------------------------------------------------------

# Load regression data
load(test_path("testdata", "resp_reg.rda"))
load(test_path("testdata", "covs_reg.rda"))

# Load classification data
load(test_path("testdata", "resp_cls.rda"))
load(test_path("testdata", "covs_cls.rda"))



# REG and CLS fitting -----------------------------------------------------
# No differences between splitting methods; use faster one for each task type

test_that("eforest works for regression", {
  expect_warning(eforest_fit <- eforest(resp_reg, covs_reg,
                                        split_type = 'coeff', ntrees = 10),
                 'No names available for covariates. Numbers are used instead.')
  expect_true(all(sapply(eforest_fit$ensemble, 
                         function(t) 'constparty' %in% class(t))))
  expect_true(all(sapply(eforest_fit$ensemble, 
                         function(t) 'party' %in% class(t))))
  expect_type(eforest_fit$oob_score, 'double')
  expect_type(eforest_fit$perf_metric, 'character')
})

test_that("eforest works for classification", {
  expect_warning(eforest_fit <- eforest(resp_cls, covs_cls, 
                                        split_type = 'cluster', ntrees = 10),
                 'No names available for covariates. Numbers are used instead.')
  expect_true(all(sapply(eforest_fit$ensemble, 
                         function(t) 'constparty' %in% class(t))))
  expect_true(all(sapply(eforest_fit$ensemble, 
                         function(t) 'party' %in% class(t))))
  expect_type(eforest_fit$oob_score, 'double')
  expect_type(eforest_fit$perf_metric, 'character')
})



# Wrong input -------------------------------------------------------------

test_that("eforest does not work with wrong input", {
  expect_error(eforest(list(resp_reg), covs_reg),
               "Argument 'response' must be provided either as a factor or as an
         object of mode 'numeric'")
  expect_error(eforest(as.character(resp_cls), covs_cls),
               "Argument 'response' must be provided either as a factor or as an
         object of mode 'numeric'")
  expect_error(eforest(resp_cls == 'Bel', covs_cls),
               "Argument 'response' must be provided either as a factor or as an
         object of mode 'numeric'")
  expect_error(eforest(resp_reg, covs_reg[[3]]),
               "Argument 'covariates' must be provided as a list")
  expect_error(eforest(resp_reg, covs_reg[[4]]),
               "Argument 'covariates' must be provided as a list")
  expect_error(eforest(resp_reg, covs_reg, minbucket = -1))
  expect_error(eforest(resp_reg, covs_reg, alpha = -3))
  expect_error(eforest(resp_reg, covs_reg, alpha = 2))
  expect_error(eforest(resp_reg, covs_reg, ntrees = -20))
  expect_error(eforest(resp_cls, covs_cls, perf_metric = 'Accuracy'),
               "The value provided for 'perf_metric' is not available for this
             task")
  expect_error(eforest(resp_reg, covs_reg, perf_metric = 'Mean Square Error'),
               "The value provided for 'perf_metric' is not available for this
             task")
})



# Right input and expected effect -----------------------------------------

test_that("input values have the expected effect", {
  
  # ntrees
  ntrees <- 5
  suppressWarnings(eforest_fit <- eforest(resp_reg, covs_reg, ntrees = ntrees))
  expect_length(eforest_fit$ensemble, ntrees)
  
  # perf_metric
  suppressWarnings(eforest_fit <- eforest(resp_reg, covs_reg, ntrees = 5,
                                          perf_metric = NULL))
  expect_identical(eforest_fit$perf_metric, 'RMSPE')
  suppressWarnings(eforest_fit <- eforest(resp_cls, covs_cls, ntrees = 5,
                                          perf_metric = NULL))
  expect_identical(eforest_fit$perf_metric, 'BAcc')
  
  # verbose
  expect_output(suppressWarnings(eforest(resp_reg, covs_reg, ntrees = 5,
                                         verbose = TRUE)))
  
})



# Predictions -------------------------------------------------------------
# No differences between splitting methods; use faster one for each task type

# Newdata for regression
new_idx <- c(76:125)
new_resp_reg <- data_reg$resp[new_idx]
new_covs_reg <- lapply(data_reg$covs, function(cov) cov[new_idx])

# Newdata for classification
new_idx <- c(51:100)
new_resp_cls <- data_cls$resp[new_idx]
new_covs_cls <- lapply(data_cls$covs, function(cov) cov[new_idx])

test_that("predict.eforest works for regression", {
  expect_warning(eforest_fit <- eforest(resp_reg, covs_reg, ntrees = 10,
                                        split_type = 'cluster'),
                 'No names available for covariates. Numbers are used instead.')
  pred <- predict(eforest_fit)
  expect_type(pred, typeof(resp_reg))
  expect_length(pred, length(resp_reg))
  pred2 <- predict(eforest_fit, newdata = new_covs_reg)
  expect_type(pred2, typeof(resp_reg))
  expect_length(pred2, length(covs_reg[[1]]))
})

test_that("predict.eforest works for classification", {
  expect_warning(eforest_fit <- eforest(resp_cls, covs_cls, ntrees = 10,
                                        split_type = 'cluster'),
                 'No names available for covariates. Numbers are used instead.')
  pred <- predict(eforest_fit)
  expect_type(pred, typeof(resp_cls))
  expect_length(pred, length(resp_cls))
  pred2 <- predict(eforest_fit, newdata = new_covs_cls)
  expect_type(pred2, typeof(resp_cls))
  expect_length(pred2, length(covs_cls[[1]]))
})




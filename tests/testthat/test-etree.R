

# Load data ---------------------------------------------------------------

# Load regression data
load(test_path("testdata", "resp_reg.rda"))
load(test_path("testdata", "covs_reg.rda"))

# Load classification data
load(test_path("testdata", "resp_cls.rda"))
load(test_path("testdata", "covs_cls.rda"))



# Test 4 + 2 combinations -------------------------------------------------

test_that("etree works for regression using 'coeff'", {
  expect_warning(etree_fit <- etree(resp_reg, covs_reg, split_type = 'coeff'),
                 'No names available for covariates. Numbers are used instead.')
  expect_s3_class(etree_fit, 'constparty')
  expect_s3_class(etree_fit, 'party')
  expect_output(print(etree_fit))
  expect_length(etree_fit$fitted$`(fitted)`, NROW(etree_fit$data[[1]]))
  expect_error(expect_output(print(etree_fit), regexp = 'n = 0'))
  #TODO: find a less convoluted way?
})


test_that("etree works for regression using 'cluster'", {
  expect_warning(etree_fit <- etree(resp_reg, covs_reg, split_type = 'cluster'),
                 'No names available for covariates. Numbers are used instead.')
  expect_s3_class(etree_fit, 'constparty')
  expect_s3_class(etree_fit, 'party')
  expect_output(print(etree_fit))
  expect_length(etree_fit$fitted$`(fitted)`, NROW(etree_fit$data[[1]]))
  expect_error(expect_output(print(etree_fit), regexp = 'n = 0'))
})


test_that("etree works for classification using 'coeff'", {
  expect_warning(etree_fit <- etree(resp_cls, covs_cls, split_type = 'coeff'),
                 'No names available for covariates. Numbers are used instead.')
  expect_s3_class(etree_fit, 'constparty')
  expect_s3_class(etree_fit, 'party')
  expect_output(print(etree_fit))
  expect_length(etree_fit$fitted$`(fitted)`, NROW(etree_fit$data[[1]]))
  expect_error(expect_output(print(etree_fit), regexp = 'n = 0'))
})


test_that("etree works for classification using 'cluster'", {
  expect_warning(etree_fit <- etree(resp_cls, covs_cls, split_type = 'cluster'),
                 'No names available for covariates. Numbers are used instead.')
  expect_s3_class(etree_fit, 'constparty')
  expect_s3_class(etree_fit, 'party')
  expect_output(print(etree_fit))
  expect_length(etree_fit$fitted$`(fitted)`, NROW(etree_fit$data[[1]]))
  expect_error(expect_output(print(etree_fit), regexp = 'n = 0'))
})


test_that("etree works for regression using 'coeff' and 'traditional", {
  expect_warning(etree_fit <- etree(resp_reg, covs_reg, split_type = 'coeff',
                                    coeff_split_type = 'traditional'),
                 'No names available for covariates. Numbers are used instead.')
  expect_s3_class(etree_fit, 'constparty')
  expect_s3_class(etree_fit, 'party')
  expect_output(print(etree_fit))
  expect_length(etree_fit$fitted$`(fitted)`, NROW(etree_fit$data[[1]]))
  expect_error(expect_output(print(etree_fit), regexp = 'n = 0'))
  #TODO: find a less convoluted way?
})


test_that("etree works for classification using 'coeff' and 'traditional", {
  expect_warning(etree_fit <- etree(resp_cls, covs_cls, split_type = 'coeff',
                                    coeff_split_type = 'traditional'),
                 'No names available for covariates. Numbers are used instead.')
  expect_s3_class(etree_fit, 'constparty')
  expect_s3_class(etree_fit, 'party')
  expect_output(print(etree_fit))
  expect_length(etree_fit$fitted$`(fitted)`, NROW(etree_fit$data[[1]]))
  expect_error(expect_output(print(etree_fit), regexp = 'n = 0'))
})



# Useful methods ----------------------------------------------------------

test_that("methods for fitted objects work", {
  
  # Fit
  expect_warning(etree_fit <- etree(resp_reg, covs_reg, split_type = 'coeff'),
                 'No names available for covariates. Numbers are used instead.')
  
  # nodeapply
  expect_output(nodeapply(etree_fit, ids = nodeids(etree_fit), print))
  
  # [[
  expect_output(expect_identical(print(etree_fit[[1]]), print(etree_fit)))
  
  # [
  expect_output(expect_identical(print(etree_fit[1]), print(etree_fit)))
  
})



# Wrong input -------------------------------------------------------------

test_that("etree does not work with wrong input", {
  expect_error(etree(list(resp_reg), covs_reg),
  "Argument 'response' must be provided either as a factor or as an
         object of mode 'numeric'")
  expect_error(etree(as.character(resp_cls), covs_cls),
  "Argument 'response' must be provided either as a factor or as an
         object of mode 'numeric'")
  expect_error(etree(resp_cls == 'Bel', covs_cls),
  "Argument 'response' must be provided either as a factor or as an
         object of mode 'numeric'")
  expect_error(etree(resp_reg, covs_reg[[3]]),
               "Argument 'covariates' must be provided as a list")
  expect_error(etree(resp_reg, covs_reg[[4]]),
               "Argument 'covariates' must be provided as a list")
  expect_error(etree(resp_reg, covs_reg, minbucket = -1))
  expect_error(etree(resp_reg, covs_reg, alpha = -3))
  expect_error(etree(resp_reg, covs_reg, alpha = 2))
})



# Right input and expected effect -----------------------------------------

test_that("minbucket and alpha have the expected effect", {
  
  # All the adjusted pvalues are less than or equal to alpha
  a <- 0.05
  suppressWarnings(etree_fit <- etree(data_reg$resp, data_reg$covs,
                                      alpha = a, p_adjust_method = "fdr"))
  pvalues <- unlist(nodeapply(etree_fit, ids = nodeids(etree_fit),
                              FUN = function(n) info_node(n)$pvalue))
  expect_true(all(p.adjust(pvalues, method = 'fdr') <= a))
  
  # All the terminal nodes have size larger than or equal to minbucket
  m <- 5
  suppressWarnings(etree_fit <- etree(data_reg$resp, data_reg$covs,
                                      minbucket = m))
  tnode_sizes <- sapply(nodeids(etree_fit, terminal = TRUE),
                        function(id) nrow(etree_fit[[id]]$data))
  expect_true(all(tnode_sizes >= m))
  
  # Minbucket equal to floor(nobs / 2) + 1 (or larger) result in a stump
  m <- floor(length(data_reg$resp) / 2) + 1
  suppressWarnings(etree_fit <- etree(data_reg$resp, data_reg$covs,
                                      minbucket = m))
  expect_length(nodeids(etree_fit), 1)
  
})

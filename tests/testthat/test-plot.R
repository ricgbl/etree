# Load data ---------------------------------------------------------------

# Load regression data
load(test_path("testdata", "resp_reg.rda"))
load(test_path("testdata", "covs_reg.rda"))

# Load classification data
load(test_path("testdata", "resp_cls.rda"))
load(test_path("testdata", "covs_cls.rda"))



# Test 4 combinations -----------------------------------------------------

test_that("plot works for regression using 'coeff'", {
  expect_warning(etree_fit <- etree(resp_reg, covs_reg, split_type = 'coeff'),
                 'No names available for covariates. Numbers are used instead.')
  expect_error(plot(etree_fit), regexp = NA)
})


test_that("plot works for regression using 'cluster'", {
  expect_warning(etree_fit <- etree(resp_reg, covs_reg, split_type = 'cluster'),
                 'No names available for covariates. Numbers are used instead.')
  expect_error(plot(etree_fit), regexp = NA)
})


test_that("plot works for classification using 'coeff'", {
  expect_warning(etree_fit <- etree(resp_cls, covs_cls, split_type = 'coeff'),
                 'No names available for covariates. Numbers are used instead.')
  expect_error(plot(etree_fit), regexp = NA)
})


test_that("plot works for classification using 'cluster'", {
  expect_warning(etree_fit <- etree(resp_cls, covs_cls, split_type = 'cluster'),
                 'No names available for covariates. Numbers are used instead.')
  expect_error(plot(etree_fit), regexp = NA)
})



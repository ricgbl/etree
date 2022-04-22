# Load data ---------------------------------------------------------------

# Load regression data
load(test_path("testdata", "resp_reg.rda"))
load(test_path("testdata", "covs_reg.rda"))
load(test_path("testdata", "new_resp_reg.rda"))
load(test_path("testdata", "new_covs_reg.rda"))

# Load classification data
load(test_path("testdata", "resp_cls.rda"))
load(test_path("testdata", "covs_cls.rda"))
load(test_path("testdata", "new_resp_cls.rda"))
load(test_path("testdata", "new_covs_cls.rda"))



# Regression (w/o newdata) ------------------------------------------------

test_that("predict works for regression using 'coeff'", {
  expect_warning(etree_fit <- etree(resp_reg, covs_reg, split_type = 'coeff'),
                 'No names available for covariates. Numbers are used instead.')
  pred <- predict(etree_fit)
  expect_type(pred, typeof(resp_reg))
  expect_length(pred, length(resp_reg))
  expect_length(unique(pred), length(nodeids(etree_fit, terminal = TRUE)))
  tnode_pred <- sapply(nodeids(etree_fit, terminal = TRUE),
                      function(id) mean(etree_fit[[id]]$`fitted`$`(response)`))
  tnode_size <- sapply(nodeids(etree_fit, terminal = TRUE),
                       function(id) sum(etree_fit[[id]]$`fitted`$`(fitted)`))
  expect_setequal(round(unique(pred), 5), round(tnode_pred, 5))
  expect_setequal(table(pred), tnode_size)
})

test_that("predict works for regression using 'cluster'", {
  expect_warning(etree_fit <- etree(resp_reg, covs_reg, split_type = 'cluster'),
                 'No names available for covariates. Numbers are used instead.')
  pred <- predict(etree_fit)
  expect_type(pred, typeof(resp_reg))
  expect_length(pred, length(resp_reg))
  expect_length(unique(pred), length(nodeids(etree_fit, terminal = TRUE)))
  tnode_pred <- sapply(nodeids(etree_fit, terminal = TRUE),
                       function(id) mean(etree_fit[[id]]$`fitted`$`(response)`))
  tnode_size <- sapply(nodeids(etree_fit, terminal = TRUE),
                       function(id) sum(etree_fit[[id]]$`fitted`$`(fitted)`))
  expect_setequal(round(unique(pred), 5), round(tnode_pred, 5))
  expect_setequal(table(pred), tnode_size)
})



# Regression (w/ data) ----------------------------------------------------

test_that("predict with newdata works for regression using 'coeff'", {
  expect_warning(etree_fit <- etree(resp_reg, covs_reg, split_type = 'coeff'),
                 'No names available for covariates. Numbers are used instead.')
  pred <- predict(etree_fit, newdata = new_covs_reg)
  expect_type(pred, typeof(resp_reg))
  expect_length(pred, length(covs_reg[[1]]))
  expect_lte(length(unique(pred)), length(nodeids(etree_fit, terminal = TRUE)))
  tnode_pred <- sapply(nodeids(etree_fit, terminal = TRUE),
                       function(id) mean(etree_fit[[id]]$`fitted`$`(response)`))
  expect_true(all(round(unique(pred), 5) %in% round(tnode_pred, 5)))
})

test_that("predict with newdata works for regression using 'cluster'", {
  expect_warning(etree_fit <- etree(resp_reg, covs_reg, split_type = 'cluster'),
                 'No names available for covariates. Numbers are used instead.')
  pred <- predict(etree_fit, newdata = new_covs_reg)
  expect_type(pred, typeof(resp_reg))
  expect_length(pred, length(covs_reg[[1]]))
  expect_lte(length(unique(pred)), length(nodeids(etree_fit, terminal = TRUE)))
  tnode_pred <- sapply(nodeids(etree_fit, terminal = TRUE),
                       function(id) mean(etree_fit[[id]]$`fitted`$`(response)`))
  expect_true(all(round(unique(pred), 5) %in% round(tnode_pred, 5)))
})



# Classification (w/o newdata) --------------------------------------------

test_that("predict works for classification using 'coeff'", {
  expect_warning(etree_fit <- etree(resp_cls, covs_cls, split_type = 'coeff'),
                 'No names available for covariates. Numbers are used instead.')
  pred <- predict(etree_fit)
  expect_type(pred, typeof(resp_cls))
  expect_length(pred, length(resp_cls))
  expect_length(unique(names(pred)), length(nodeids(etree_fit, terminal = TRUE)))
  tnode_pred <- factor(sapply(nodeids(etree_fit, terminal = TRUE), function(id)
    names(which.max(table(etree_fit[[id]]$`fitted`$`(response)`)))))
  tnode_size <- sapply(nodeids(etree_fit, terminal = TRUE), function(id)
    sum(etree_fit[[id]]$`fitted`$`(fitted)`))
  expect_setequal(pred[unique(names(pred))], tnode_pred)
  expect_setequal(table(names(pred)), tnode_size)
})

test_that("predict works for classification using 'cluster'", {
  expect_warning(etree_fit <- etree(resp_cls, covs_cls, split_type = 'cluster'),
                 'No names available for covariates. Numbers are used instead.')
  pred <- predict(etree_fit)
  expect_type(pred, typeof(resp_cls))
  expect_length(pred, length(resp_cls))
  expect_length(unique(names(pred)), length(nodeids(etree_fit, terminal = TRUE)))
  tnode_pred <- factor(sapply(nodeids(etree_fit, terminal = TRUE), function(id)
    names(which.max(table(etree_fit[[id]]$`fitted`$`(response)`)))))
  tnode_size <- sapply(nodeids(etree_fit, terminal = TRUE), function(id)
    sum(etree_fit[[id]]$`fitted`$`(fitted)`))
  expect_setequal(pred[unique(names(pred))], tnode_pred)
  expect_setequal(table(names(pred)), tnode_size)
})



# Classification (w/ data) ------------------------------------------------

test_that("predict with newdata works for classification using 'coeff'", {
  expect_warning(etree_fit <- etree(resp_cls, covs_cls, split_type = 'coeff'),
                 'No names available for covariates. Numbers are used instead.')
  pred <- predict(etree_fit, newdata = new_covs_cls)
  expect_type(pred, typeof(resp_cls))
  expect_length(pred, length(covs_reg[[1]]))
  expect_lte(length(unique(names(pred))), length(nodeids(etree_fit, terminal = TRUE)))
  tnode_pred <- factor(sapply(nodeids(etree_fit, terminal = TRUE), function(id)
    names(which.max(table(etree_fit[[id]]$`fitted`$`(response)`)))))
  expect_true(all(pred[unique(names(pred))] %in% tnode_pred))
})

test_that("predict with newdata works for classification using 'cluster'", {
  expect_warning(etree_fit <- etree(resp_cls, covs_cls, split_type = 'cluster'),
                 'No names available for covariates. Numbers are used instead.')
  pred <- predict(etree_fit, newdata = new_covs_cls)
  expect_type(pred, typeof(resp_cls))
  expect_length(pred, length(covs_reg[[1]]))
  expect_length(unique(names(pred)), length(nodeids(etree_fit, terminal = TRUE)))
  tnode_pred <- factor(sapply(nodeids(etree_fit, terminal = TRUE), function(id)
    names(which.max(table(etree_fit[[id]]$`fitted`$`(response)`)))))
  expect_true(all(pred[unique(names(pred))] %in% tnode_pred))
})


test_that("dist_comp works well with logical vectors", {
  nobs <- 10
  obj <- as.logical(rbinom(nobs, 1, 0.5))
  expect_error(d <- dist_comp(obj), regexp = NA)
  expect_s3_class(d, 'dist')
  expect_length(d, nobs * (nobs - 1) / 2)
})

test_that("dist_comp works well with integer vectors", {
  nobs <- 10
  obj <- rpois(nobs, 5)
  expect_error(d <- dist_comp(obj), regexp = NA)
  expect_s3_class(d, 'dist')
  expect_length(d, nobs * (nobs - 1) / 2)
})

test_that("dist_comp works well with numeric vectors", {
  nobs <- 10
  obj <- rnorm(nobs)
  expect_error(d <- dist_comp(obj), regexp = NA)
  expect_s3_class(d, 'dist')
  expect_length(d, nobs * (nobs - 1) / 2)
})

test_that("dist_comp works well with factors", {
  nobs <- 10
  obj <- factor(letters[1:nobs])
  expect_error(d <- dist_comp(obj), regexp = NA)
  expect_s3_class(d, 'dist')
  expect_length(d, nobs * (nobs - 1) / 2)
})

test_that("dist_comp works well with functional data", {
  nobs <- 10
  obj <- fda.usc::rproc2fdata(nobs, seq(0, 1, len = 100), sigma = 1)
  expect_error(d <- dist_comp(obj), regexp = NA)
  expect_s3_class(d, 'dist')
  expect_length(d, nobs * (nobs - 1) / 2)
})

test_that("dist_comp works well with graphs", {
  nobs <- 10
  obj <- lapply(1:nobs, function(j) igraph::sample_gnp(100, 0.2))
  expect_error(d <- dist_comp(obj), regexp = NA)
  expect_s3_class(d, 'dist')
  expect_length(d, nobs * (nobs - 1) / 2)
})

test_that("dist_comp works well with directed graphs", {
  nobs <- 10
  obj <- lapply(1:nobs, function(j) igraph::sample_gnp(100, 0.2, directed = T))
  expect_error(d <- dist_comp(obj), regexp = NA)
  expect_s3_class(d, 'dist')
  expect_length(d, nobs * (nobs - 1) / 2)
})

test_that("dist_comp works well with directed graphs", {
  nobs <- 10
  obj <- lapply(1:nobs, function(j) {
    g <- igraph::sample_gnp(100, 0.2, directed = T)
    igraph::E(g)$weight <- runif(length(igraph::E(g)), 0, 1)
    return(g)
    })
  expect_error(d <- dist_comp(obj), regexp = NA)
  expect_s3_class(d, 'dist')
  expect_length(d, nobs * (nobs - 1) / 2)
})

test_that("dist_comp works well with persistence diagrams", {
  nobs <- 10
  x <- lapply(rep(100, nobs), function(np) TDA::circleUnif(np))
  obj <- lapply(x, TDA::ripsDiag, maxdimension = 1, maxscale = 3)
  expect_error(d <- dist_comp(obj), regexp = NA)
  expect_s3_class(d, 'dist')
  expect_length(d, nobs * (nobs - 1) / 2)
})

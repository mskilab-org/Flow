## tests for Flow.R

library(Flow)

context('test Flow operations')


test_that("module", {
    expect_equal(2+2, 4)
})
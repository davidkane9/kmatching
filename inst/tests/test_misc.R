context("Miscelaneous Tests")

test_that("Dummy works", {
  
  ## Not sure why these should produce different matrices . . .
  
  expect_that(dummy(letters[1:3]), equals(matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3)))
  expect_that(dummy(letters[3:1]), equals(matrix(c(0,0,1,0,1,0,1,0,0), nrow = 3)))
})



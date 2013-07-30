context("Miscelaneous Tests")

test_that("Dummy works", {
  
  ## Note that dummy sorts the incoming vector.
  
  expect_that(dummy(letters[1:3]), equals(matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3)))
  expect_that(dummy(letters[3:1]), equals(matrix(c(0,0,1,0,1,0,1,0,0), nrow = 3)))
  
  ## Not sure what the appropriate NA behavior is, but wanted to test current
  ## output.
  
  expect_that(dummy(c(letters[1:3], NA)), equals(matrix(c(1,0,0,0,1,0,0,0,1,0,0,0), nrow = 3)))
  
  
})



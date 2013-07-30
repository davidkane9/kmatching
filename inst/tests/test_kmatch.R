context("Test kmatch")

## Ought to consider pulling out Lalonde data and adding it to the kmatching
## package, but not today.


test_that("Using Lalonde produces results that match perfectly", {
  library(MatchIt)
  data(lalonde)
  
  ## Single 0/1 variable.
  
  z <- kmatch(lalonde, match.var = c("hispan"), weight.var = "treat", n = 10)
    
  for(i in 1:10){
  expect_that(sum(lalonde$hispan[lalonde$treat == 1]), 
              equals(apply(lalonde$hispan * z, 2, sum)[i]))
  }  
  
  ## Double check that it works for a factor variable as well.
  
  lalonde2 <- lalonde
  lalonde2$hispan <- factor(lalonde2$hispan, labels = c("Non-hispanic", "Hispanic"))
  
  z <- kmatch(lalonde2, match.var = c("hispan"), weight.var = "treat", n = 10)
  for(i in 1:10){
    expect_that(sum(lalonde$hispan[lalonde$treat == 1]), 
                equals(apply(lalonde$hispan * z, 2, sum)[i]))
  } 
  
  ## Two 0/1 variables.
  
  z <- kmatch(lalonde, match.var = c("hispan", "black"), weight.var = "treat", n = 10)
  
  for(i in 1:10){
    expect_that(sum(lalonde$hispan[lalonde$treat == 1]), 
                equals(apply(lalonde$hispan * z, 2, sum)[i]))
    expect_that(sum(lalonde$black[lalonde$treat == 1]), 
                equals(apply(lalonde$black * z, 2, sum)[i]))
  }
   
  ## Two 0/1 variables and continuous variable.
  
  z <- kmatch(lalonde, match.var = c("hispan", "black", "educ"), weight.var = "treat", n = 10)
  
  ## Be careful in the use of sum versus mean. You want 'mean' to work as well,
  ## but it doesn't because the number of observation is different when
  ## comparing the subset which are Hispanic to a matching vector which includes
  ## the entire sample.
  
  ## Side note. This is an interesting example because it highlights the 
  ## slow-mixing nature of the process. The 10 solution vectors are virtually 
  ## identical. Even if n = 100, it is hard to tell the samples apart. You need 
  ## n = 1,000 or, better, 10,000 to start to see meaningful differences in the 
  ## samples. Is there a way to speed up that process? How can we know if the 
  ## mixing is complete? I think that rstan might provide some guidance here. We
  ## want kmatch to allow for --- even to do so by default --- starting 4 or
  ## more different chains, from dispersed locations, and then combining the
  ## results.
  
  ## Another interesting thing is that, even when I run the above command for n 
  ## = 100,000, all the samples have the same number of non-zero weights: 429. I
  ## guess that this makes sense in that there are 429 control subjects to 
  ## choose from, but there are obviously answers that have non-zero weights on 
  ## fewer than 429 of them. Indeed, it seems like anywhere from 150 to 225 have
  ## very low, almost zero, weight placed on them. Might be nice to encourage
  ## the process to just zero those units out and then reweight the rest, if
  ## only for aesthetic reasons.
  
  ## It is also interesting to look at specific individuals. Considers these two:
  
  ## > lalonde[200:201,]
  ## treat age educ black hispan married nodegree      re74     re75      re78
  ## PSID15     0  22   14     1      0       1        0  748.4399 11105.37 18208.550
  ## PSID16     0  42    0     0      1       1        1 2797.8330 10929.92  9922.934
  
  ## Because the average education for the treated is 10 and the average age is 
  ## 26, PSID15 is a good match. In all 100,000 samples, his weighting ranges 
  ## from 1.6 to 2.7. PSID16, on the other hand, is very different from the
  ## center of mass of the treated, so his weighting is mostly below 0.1 and
  ## never above 0.25. 
  
  ## Another fun activity is to look at the 100,000 sample weights for 
  ## individual units like PSID15 and PSID16. Looking at a histogram can show, 
  ## among (I assume) many possible patterns, a peak near zero with an 
  ## exponential drop of to toward higher weights (PSID16) or a bi-modal 
  ## structure (PSID15). It is also interesting to consider the time series of 
  ## weights. We want these to show mixing. Often they do, especially with 
  ## skiplength = 10 or 100. But not always! Even with skiplength = 1,000, we
  ## don't get even reasonable mixing for row 200 with matchvars  = c("age",
  ## "educ", "black").
  
  
  for(i in 1:10){
    expect_that(sum(lalonde$hispan[lalonde$treat == 1]), 
                equals(apply(lalonde$hispan * z, 2, sum)[i]))
    expect_that(sum(lalonde$black[lalonde$treat == 1]), 
                equals(apply(lalonde$black * z, 2, sum)[i]))
    expect_that(sum(lalonde$educ[lalonde$treat == 1]), 
                equals(apply(lalonde$educ * z, 2, sum)[i]))
  }
  
})

test_that("Mixing properly", {
  matchvars  = c("age", "educ", "black")
  k = kmatch(x = lalonde, weight.var  = "treat", match.var = matchvars, n = 1000)
  
})

## Add testcase with NAs



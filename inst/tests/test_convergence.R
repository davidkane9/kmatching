
context("Testing for convergence with Coda")

## fake data
treated = rep(0, 500)
treated[1:100] = 1
gender = rep("male", 500)
gender[sample(1:500, 250)] = "female"
bodyfat.zscore = rnorm(500)
smokes = rep("no", 500)
## high smoking treatment group
smokes[sample(1:100, 40)] = "yes"
smokes[sample(101:500, 80)] = "yes"


heart.disease = data.frame(treated = treated, gender = gender, bodyfat = bodyfat.zscore, smokes = smokes)

matchvars = c("gender", "bodyfat", "smokes")
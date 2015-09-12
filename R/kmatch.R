#' Create matching vectors based on data 
#' 
#' Given a data frame and the specified matching vectors and weight vector, 
#' \code{kmatch} extracts the matching vectors and weight vector from the 
#' data set, and form them to a system of linear equations, whose matrix
#' expression is \eqn{Aw=b}, where A is the matching vectors and b is the
#' weight vector. The \eqn{w} is the weights we want. 
#' \code{kmatch} then randomly samples \eqn{w} subject to the
#' non-negativity constraint \eqn{w_i >= 0} and \eqn{sum(w)=0}. These two
#' make the sampling space bounded so
#' uniform sampling over the whole sampling space is achievable. Then the
#' \eqn{w} is the weights we want.
#' 
#' @param x is a data frame containing needed input data
#' @param weight.var character name of the column of the input weights
#' @param match.var character vector of names of columns of 'data' we wish to match on
#' @param n numeric number of weight vectors desired. Default is 1
#' @param thin only 1 out of every \code{thin} vectors will be recorded. Default is 5
#' @param replace if FALSE then any rows that had non-zero weight in \code{weight.var}
#'        will be set to 0 in the output. Default is \code{FALSE}
#' @param chains number of chains. Default is 1.
#' @param burn is the number of points to burn in hitandrun (not very costly, because it just adds to total number
#'        of iterations, whereas thin is a multiple)
#' 
#' @return Returns a list of "chains" matrices of 'n' sets of weights that match the given set of
#' weights in terms of weighted averages to the 'match.var' factors. The columns of each matrix
#' are the sets of weights.
#' @export
#' @author Mike Flynn \email{mflynn210@@gmail.com}
#' @examples
#' set.seed(40)
#' x <- data.frame(size = rnorm(10), weight = rep(.02, 10))
#' weights <- kmatch(x, weight.var = "weight", match.var = "size", n = 100, replace = TRUE)
#' 
#' x2 <- data.frame(var1 = c(-2,-1, 0, 1, 2), 
#'                  var2 = c(1, -1, 2, -1, 0), 
#'                  w = rep(.2, 5))
#' weights2 <- kmatch(x = x2, weight.var = "w", 
#'                    match.var = c("var1", "var2"), 
#'                    n = 10, thin = 1, replace= TRUE)[[1]]

kmatch <- function(x,
                   weight.var,
                   match.var,
                   n       = 1,
                   thin    = 5,
                   replace = FALSE,
                   chains  = 1,
                   burn    = 0.5) {
  
  ## Input checking.
  
  stopifnot(is.data.frame(x))
  stopifnot(is.character(weight.var))
  stopifnot(is.character(match.var))
  stopifnot(n > 0)
  stopifnot(all(c(weight.var, match.var) %in% names(x)))
  
  notincluded <- numeric()
  ret <- matrix(0, nrow = nrow(x), ncol = n)
  
  if(any(is.na(x[[weight.var]]))) {
    warning(paste(length(which(is.na(x[[weight.var]]))), " entries of '", 
                  weight.var, "' are missing, these rows are zeroed in the output", 
                  sep = ""))
    notincluded <- c(notincluded, which(is.na(x[[weight.var]])))
    x <- subset(x, !is.na(x[[weight.var]]))
  }
  
  for(v in match.var) {
    if(any(is.na(x[[v]]))) {
      warning(paste(length(which(is.na(x[[v]]))), " entries of '", v, "' 
                    are missing, these rows are zeroed in the output", sep = ""))
      notincluded <- c(notincluded, which(is.na(x[[v]])))
      x <- subset(x, !is.na(x[[v]]))
    }
  }
  
  ## First step is to build the constraint equation that samples will need to 
  ## satisfy: Aw = b. Note that this equation does not include the inequality 
  ## constraint that all elements of x must be >= 0. A is a matrix which 
  ## includes the match.var variable exposures as rows. Multiplying these rows 
  ## by the column vector x equals the appropriate element in b. 
  
  ## The number of constraints, and hence the number of rows in A and elements 
  ## in b, equals the number of levels of all factor/character variables in 
  ## match.var plus the number of numeric variables plus one. The last row is 
  ## due to the requirement that the sum of the output weights must equal the 
  ## sum of the input weight.var. So, if match.var includes country (factor 
  ## variable with 10 levels), size and value (both numeric), there would be 13
  ## rows in A and elements in b.
  
  equation <- constraint_equation(x, weight.var, match.var, replace)
  
  ## Generate the weights. Are we dealing with number of chains in a sensible
  ## way? Surely, this needs to be a choice variable for kmatch.
  #print(ncol(equation$A))
  
  ## don't need to care abt linear dependency because walkr takes care of that 

  
  if(is.matrix(equation$A)) {
    print(equation$A)
    print(equation$b)
    weights <- walkr::walkr(A = equation$A, b = equation$b, 
                            points = n, thin = thin, 
                            chains = chains, burn = burn)
  }
  
  ## handle the case where A only has 1 row, so then R automatically 
  ## converts it to a numeric
  
  else {
    print(matrix(equation$A, ncol = length(equation$A)))
    print(equation$b)
    return(list(matrix(equation$A, ncol = length(equation$A)), equation$b))
    weights <- walkr::walkr(A = matrix(equation$A, ncol = length(equation$A)), b = equation$b, 
                            points = n, thin = thin, 
                            chains = chains, burn = burn)
    
  }
  
  
  #return(weights)
    
    
    ## old code
    ## hitandrun(A = equation$A, b = equation$b, n = n, thin = thin, 
    ##                    check = check, chains = chains, burn = burn)
  
  ## Now that we have the weights, we create the return matrix, taking account of
  ## the value of replace.
  
  m <- matrix(0, ncol = ncol(weights), nrow = nrow(x))     
  
  if(!replace) {
    m[which(x[[weight.var]] == 0),] <- weights
  } else {
    m <- weights
  }
  
  if(length(notincluded) > 0) {
    ret[-notincluded,] <- m
  } else {
    ret <- m
  }
  
  weights <- ret
  
  return(weights)
}

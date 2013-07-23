#' Create matching vectors based on data Creates a matrix of weights which
#' 
#' @param x data frame containing needed input data
#' @param weight.var character name of the column of the input weights
#' @param match.var character vector of names of columns of 'data' we wish to
#'   match on
#' @param n numeric number of weight vectors desired. Default is 1
#' @param replace logical indicating whether or not bservations weighted in the
#'   original weight.var are allowed positive weight in the output. Default is
#'   FALSE
#' @param ... parameters to be passed to the sampling methods
#' @export
#' @author Mike Flynn \email{mflynn210@@gmail.com}
#'   
#'   
#' @examples
#' x <- data.frame(size = rnorm(50), weight = rep(.02, 50))
#' weights <- kmatch(x, weight.var = "weight", match.var = "size", n = 100, replace = TRUE)

kmatch <- function(x, weight.var, match.var,  n = 1, replace = FALSE, ...) {
    
  stopifnot(n > 0)
  stopifnot(all(c(weight.var, match.var) %in% names(x)))
  
  ## Intialize list that will be turned into a matrix with do.call
  
  Alist = list()
    
  ## Include the continuous and discrete variables in Alist. Need to be careful
  ## in deciding just what a 0/1 variable means, for example.
  
  for(i in 1:length(match.var)){
    if(is.numeric(x[[match.var[i]]])){
      ## continuous
      Alist[[i]] = x[[match.var[i]]]

    } else {
      ## discrete
      Alist[[i]] = dummy(x[[match.var[i]]])
    }
  }
  ## do.call on Alist makes constraint matrix A

  A = do.call(rbind, Alist)

  ## b is the constraint matrix
  b = A %*% x[[weight.var]]
  
  ## attach the "match sum" constraint, redundanies no longer matter
  sumlimit = sum(x[[weight.var]])
  A = rbind(A, rep(1, ncol(A)))
  b = c(b, sumlimit)
  
  if(!replace) {
    if(sum(x[[weight.var]] > 0) >= nrow(x) ) stop("All rows are weighted, set replace = TRUE.")
    ## remove columns corresponding to variables that have weight in the original 
    A = A[,-which(x[[weight.var]] > 0)]
  }
  
  
  weights = hitandrun(A, b, n=n, ...)
  
  ret = matrix(0, nrow = nrow(x), ncol = n)
  if(!replace) {
    ret[which(x[[weight.var]] == 0),] = weights
  } else {
    ret = weights
  }  
  return(ret)
}

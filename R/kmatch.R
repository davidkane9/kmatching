#' Create matching vectors based on data
#' 
#' An abstraction on top of the random sampling methods. This function
#' allows a user to specify which variables they would like to match exposures
#' to an original vector of weights.
#' 
#' @param data A data frame containing data we want to match
#' @param match.var A list of names of columns of 'data' we wish to match
#' @param weight.var The name of the column of orignal set of weights.
#' @param n The number of outputs desired 
#' @param ... Parameters to be passed to the sampling methods
#' @export
#' @author Mike Flynn \email{mflynn210@@gmail.com}
#' 
#' @examples
#' data = data.frame(size = rnorm(50), weight = rep(.02, 50))
#' weights = kmatch(data = data, match.var = "size", weight.var = "weight", n = 100)
#' 
kmatch <- function(data, match.var, weight.var, n, ...) {
  
  ## Where are the comments? This code is very hard to follow.
  
  Alist = list()
  hasdiscrete = FALSE
  for(i in 1:length(match.var)) {
    if(class(data[[match.var[i]]]) == "numeric") {
      Alist[[i]] = data[[match.var[i]]]
    } else {
      hasdiscrete = TRUE
      Alist[[i]] = .dummy(data[[match.var[i]]])
    }
  }
  A = do.call(rbind, Alist)
  b = A %*% data[[weight.var]]
  
  if(!hasdiscrete) {
    sumlimit = sum(data[[weight.var]])
    A = rbind(A, rep(1, ncol(A)))
    b = c(b, sumlimit)
  }
  
  weights = hitandrun(A, b, n=n, ...)
  return(weights)
}

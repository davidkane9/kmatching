#' Create matching vectors based on data
#' 
#' An abstraction on top of the random sampling methods. This function
#' allows a user to specify which variables they would like to match exposures
#' to an original vector of weights.
#' 
#' @param data A data frame containing data we want to match
#' @param match.var A list of names of columns of 'data' we wish to match
#' @param weight.var The name of the column of orignal set of weights.
#' @param sumlimit An optional parameter specifying what the weights in the vector
#' should sum to. This cannot be used if there are discrete matching variables because
#' it ends up being redundant or contradictory. 
#' @param ... Parameters to be passed to the sampling methods
#' @export
#' @author Mike Flynn \email{mjf2@@williams.edu}
#' 
#' @examples
#' data = data.frame(size = rnorm(50), weight = rep(.02, 50))
#' weights = kmatch(data = data, match.var = "size", weight.var = "weight", n = 100)
#' 
kmatch <-
function(data, match.var, weight.var, sumlimit = NULL, ...) {
  
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
  
  if(!is.null(sumlimit)) {
    if(hasdiscrete) stop("Setting sumlimit is redundant with a discrete matching variable")
    A = rbind(A, rep(1, ncol(A)))
    b = c(b, sumlimit)
  }
  weights = hitandrun(A, b, ...)
  return(weights)
}

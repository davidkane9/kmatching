#' Create matching vectors based on data
#' 
#' An abstraction on top of the random sampling methods. This function
#' allows a user to specify which variables they would like to match exposures
#' to an original vector of weights. Returns a matching portfolio the same size as
#' the original.
#' 
#' @param data A data frame containing data we want to match
#' @param match.var A list of names of columns of 'data' we wish to match
#' @param weight.var The name of the column of orignal set of weights.
#' @param n The number of outputs desired 
#' @param replace FALSE if nothing weighted in the original weight.var
#' should be weighted in the outputs.
#' @param ... Parameters to be passed to the sampling methods
#' @export
#' @author Mike Flynn \email{mflynn210@@gmail.com}
#' 
#' 
#' @examples
#' data = data.frame(size = rnorm(50), weight = rep(.02, 50))
#' weights = kmatch(data = data, match.var = "size", weight.var = "weight", n = 100, replace = TRUE)
#' 
kmatch <- function(data, match.var, weight.var, n, replace = FALSE, ...) {
  ## index will help us keep track of subsetted data
  data$index = 1:nrow(data)
  
  ## Intialize list that will be turned into a matrix with do.call
  Alist = list()
  
  ## OLD CODE FROM TRYING TO SHRINK, IGNORE #####################
#   ## seperate into discrete and continuous variables
#   match.var.c = match.var[sapply(match.var, function(x) class(data[[x]]) == "numeric" )]
#   match.var.d = match.var[sapply(match.var, function(x) class(data[[x]]) != "numeric" )]
#   
#   if(length(match.var.d) > 0) {
#     ## get a small subset of the data frame that preserves number of rows
#     ## in each variable category.
#     newdata = ddply(data, match.var.d, function(x) {
#       pos = 1:nrow(x)
#       if(!replace) pos = pos[-which(x[[weight.var]] > 0)]
#       if(length(pos) < length(which(x[[weight.var]] > 0))) {
#         stop("Error: there are not enough entries in each possible
#            cross-section of discrete variables to not weight at least one of the rows in the 
#            original data. Try setting replace = TRUE.")
#       }
#       len = length(which(x[[weight.var]] != 0))
#       return(x[sample(pos, len),])
#     })
#   } else {
#     pos = 1:nrow(data)
#     if(!replace) pos = pos[-which(data[[weight.var]] > 0)]
#     newdata = data[sample(pos, length(which(data[[weight.var]] > 0)))]
#   }
  ###########################################################
  
  ## checking if there is a discrete variable will prevent redundancies
  
  Alist2 = list()
  ## include the continuous and discrete variables in Alist
  for(i in 1:length(match.var)) {
    if(class(data[[match.var[i]]]) == "numeric") {
      ## continuous
      Alist[[i]] = data[[match.var[i]]]
      ##Alist2[[i]] = newdata[[match.var[i]]]
    } else {
      ## discrete
      Alist[[i]] = .dummy(data[[match.var[i]]])
      ##Alist2[[i]] = .dummy(newdata[[match.var[i]]])
    }
  }
  ## do.call on Alist makes constraint matrix A
  A = do.call(rbind, Alist)
  ##A2 = do.call(rbind, Alist2)
  ## b is the constraint matrix
  b = A %*% data[[weight.var]]
  
  ## attach the "match sum" constraint, redundanies no longer matter
  sumlimit = sum(data[[weight.var]])
  A = rbind(A, rep(1, ncol(A)))
  b = c(b, sumlimit)
  
  if(!replace) {
    if(sum(data[[weight.var]] > 0) >= nrow(data) ) stop("All rows are weighted, set replace = TRUE.")
    ## remove columns corresponding to variables that have weight in the original 
    A = A[,-which(data[[weight.var]] > 0)]
  }
  
  weights = hitandrun(A, b, n=n, ...)
  
  ret = matrix(0, nrow = nrow(data), ncol = n)
  if(!replace) {
    ret[which(data[[weight.var]]) == 0,] = weights
  } else {
    ret = weights
  }  
  return(weights)
}

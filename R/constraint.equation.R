#' Creates a constraint equation based on the input data frame.
#' 
#' Given an input data frame, extract the matching vectors and the weight vectors out
#' of the data frame. Then format the extracted data to form a linear constraint \eqn{Aw=b}
#' that defines a sample space from which we will later sample from.
#' 
#' @param x data frame containing needed input data
#' @param weight.var character name of the column of the input weights
#' @param match.var character vector of names of columns of 'data' we wish to 
#'        match on
#' @param replace logical indicating whether or not observations weighted in the 
#'        original weight.var are allowed positive weight in the output.
#'   
#' @return A list with two named components: A and b, representing the
#'         components of the constraint equation \eqn{Aw = b}
#'   
#' @author David Kane \email{<dave.kane@@gmail.com>}
#' @export

constraint_equation <- function(x, weight.var, match.var, replace){
  
  ## Intialize list that will be turned into a matrix with do.call
  
  Alist <- list()
  
  ## Include the continuous and discrete variables in Alist. Need to be careful
  ## in deciding just what a 0/1 variable means, for example.
  
  for(i in 1:length(match.var)){
    
    if(is.numeric(x[[match.var[i]]])){
      
      Alist[[i]] <- x[[match.var[i]]]
      
    } else {
      
      Alist[[i]] <- dummy(x[[match.var[i]]])
    }
  }
  
  ## Constructs constraint matrix (left-hand-side) A
  
  A <- do.call(rbind, Alist)
  
  ## Construct constraint matrix (right-hand-side) b
  b <- A %*% x[[weight.var]]
  
  ## attach the "match sum" constraint, redundanies no longer matter
  
  sumlimit <- sum(x[[weight.var]])
  A <- rbind(A, rep(1, ncol(A)))
  b <- c(b, sumlimit)
  
  if(!replace) {
    
    if(sum(x[[weight.var]] > 0) >= nrow(x)){
      stop("All rows are weighted, set replace = TRUE.")
    }
    
    ## remove columns corresponding to variables that have weight in the original 
    
    A <- A[, -which(x[[weight.var]] > 0)]
  }
  
  ## DO NOT USE ARROWS HERE! You are passing arguments to the list function!
  ## If you change that to arrows, there'll be notorious bugs!
  
  ## just a little hack for now to get rid of sum(x) != 1 row 
  ## should probably fix later to be better
  
  A <- A[1:(nrow(A)-1),     ]
  b <- b[1:(length(b)-1)]
  
  return(list(A = A, b = b)) 
}
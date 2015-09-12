#' A package to generate weights according to matching vectors
#' 
#' \tabular{ll}{
#' Package \tab kmatching\cr
#' Type: \tab Package\cr
#' Version: \tab 0.3-2\cr
#' Date: \tab 2015-08-28\cr
#' License: \tab GPL-3\cr
#' }
#' 
#' Given a data frame and the specified matching vectors and weight vector, 
#' the package can extract the matching vectors and weight vector from the 
#' data set, and form them to a system of linear equations, whose matrix
#' expression is \eqn{Aw=b}, where A is the matching vectors and b is the
#' weight vector. The package then randomly sample \eqn{w} subject to the
#' non-negativity constraint \eqn{w_i >= 0} and \eqn{sum(w)=0}. These two
#' constrains bound the simplex and make the sampling space bounded so
#' uniform sampling over the whole sampling space is achievable. Then the
#' \eqn{w} the package generates is the weights.
#' 
#' @name kmatching-package
#' @docType package
NULL
#> NULL
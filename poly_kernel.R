poly_kernel <- function(A, B, d) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  n <- nrow(A)
  p <- ncol(B)

  if(ncol(A) == nrow(B)){
    K <- (A%*%B)^d
  }else{
    return("ERROR: non-conformable arguments")
  }
  return(K)
}

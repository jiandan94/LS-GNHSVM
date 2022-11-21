lsgnhsvm <- function(trainx, trainy, lam1, lam2, kernel = "linear", ...){
  require("MASS")
  # compute the kernel
  if(kernel == "linear"){
    x_pos <- trainx[which(trainy == 1), ]
    x_neg <- trainx[which(trainy == -1), ]
  }else if(kernel == "gaussian"){
    xk_pos <- trainx[which(trainy == 1), ]
    xk_neg <- trainx[which(trainy == -1), ]
    xk <- rbind(xk_pos, xk_neg)
    x_pos <- gaussian_kernel(xk_pos, t(xk), ...)
    x_neg <- gaussian_kernel(xk_neg, t(xk), ...)
  }else if(kernel == "poly"){
    xk_pos <- trainx[which(trainy == 1), ]
    xk_neg <- trainx[which(trainy == -1), ]
    xk <- rbind(xk_pos, xk_neg)
    x_pos <- poly_kernel(xk_pos, t(xk), ...)
    x_neg <- poly_kernel(xk_neg, t(xk), ...)
  }
  
  n_pos <- nrow(x_pos)
  n_neg <- nrow(x_neg)
  n <- n_pos + n_neg
  p <- ncol(x_pos)
  
  e_pos <- matrix(1, nrow = n_pos, ncol = 1)
  e_neg <- matrix(1, nrow = n_neg, ncol = 1)
  e <- matrix(1, nrow = n, ncol = 1)
  
  x_tilde_pos <- cbind(x_pos, e_pos)
  x_tilde_neg <- cbind(x_neg, e_neg)
  x_tilde <- rbind(x_tilde_pos, -x_tilde_neg)
  
  # compute matrix Q
  E <- diag(p + 1)
  mat1 <- ginv(E + lam1*t(x_tilde_pos)%*%x_tilde_pos)
  mat2 <- ginv(E + lam1*t(x_tilde_neg)%*%x_tilde_neg)
  Q <- x_tilde%*%(mat1 + mat2)%*%t(x_tilde)
  
  # compute u 
  u <- ginv(Q + diag(n)/lam2)%*%e
  
  # compute w+, b+ and w-, b-
  w_tilde_pos <- mat1%*%t(x_tilde)%*%u
  w_tilde_neg <- mat2%*%t(x_tilde)%*%u
  w_pos <- w_tilde_pos[1:p]
  b_pos <- w_tilde_pos[p+1]
  w_neg <- w_tilde_neg[1:p]
  b_neg <- w_tilde_neg[p+1]
  
  gnhClassifier <- list(wpos = w_pos, bpos = b_pos, wneg = w_neg, bneg = b_neg, u = u)
  return(gnhClassifier)
}

varimax = function (x, normalize = TRUE, eps = 1e-05) 
{
  nc <- ncol(x)
  if (nc < 2) 
    return(x)
  if (normalize) {
    sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
    x <- x/sc
  }
  p <- nrow(x)
  TT <- diag(nc)
  d <- 0
  for (i in 1L:1000L) {
    z <- x %*% TT
    B <- t(x) %*% (z^3 - z %*% diag(drop(rep(1, p) %*% z^2))/p)
    sB <- La.svd(B)
    TT <- sB$u %*% sB$vt
    dpast <- d
    d <- sum(sB$d)
    if (d < dpast * (1 + eps)) 
      break
  }
  z <- x %*% TT
  if (normalize) 
    z <- z * sc
  dimnames(z) <- dimnames(x)
  class(z) <- "loadings"
  list(loadings = z, rotmat = TT)
}

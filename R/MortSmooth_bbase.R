MortSmooth_tpower <- function (x, t, p){
  (x - t)^p * (x > t)
}
MortSmooth_bbase <- function (x, xl, xr, ndx, deg){
  dx <- (xr - xl)/ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, MortSmooth_tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1)/(gamma(deg + 1) * dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  B
}
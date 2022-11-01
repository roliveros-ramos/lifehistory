
#' @export
QBD = function(x, ...) {
  UseMethod("QBD")
}

#' @export
QBD.elefan = function(x, a, b, a50, a75, age0, parallel=TRUE, control=list(), ...) {


  out = .QBD_fit(x=x$bycohort,
                 lifespan=x$lifespan, a=a, b=b, a50=a50, a75=a75, age0=age0,
                 parallel=parallel, control=control)

  opt = .QBD_par(Linf=x$par$Linf, k=x$par$k, t0=x$par$t0, a=a, b=b, beta=out$beta,
                 lifespan=x$lifespan, a50=a50, a75=a75, age0=age0)

  x$bycohort = out$bycohort
  x$a = a
  x$b = b
  x$a50 = a50
  x$a75 = a75
  x$beta = out$beta

  x$par$c = as.numeric(opt[1])
  x$par$r = as.numeric(opt[2])
  x$par$beta = out$beta

  return(x)

}

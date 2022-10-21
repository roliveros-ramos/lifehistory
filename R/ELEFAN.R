
#' @export
ELEFAN1 = function(x, ...) {
  UseMethod("ELEFAN1")
}

#' @export
ELEFAN1.len_freq = function(x, Linf, k, lifespan, birth="01-apr", method="rsa",
                            parallel=TRUE, control=list()) {

  marks = x$marks
  bins  = x$bins
  date  = x$date
  x     = x$x

  out = .ELEFAN1(x=x, marks=marks, bins=bins, date=date, Linf=Linf, k=k,
                 lifespan=lifespan, birth=birth, method=method,
                 parallel=parallel, control=control)

  class(out) = "elefan"

  return(out)

}

.ELEFAN1 = function(x, marks, bins, date, Linf, k, lifespan, birth,
                    method, parallel, control) {

  method = match.arg(method, choices = c("rsa"))

  .t0 = Sys.time()

  control = .elefan_control_defaults(control)

  n = control$n

  Linf = .check_Linf(Linf)
  k    = .check_k(k)

  # ---- OBSERVED
  rcounts  = lf_restructure(x, MA=control$MA, method=control$method)
  peaks    = ASP(rcounts, return.peaks = TRUE)
  # new form of observed data
  obs = expand.grid(length=marks, date=date)
  obs$Fs = as.numeric(rcounts)  # turn matrix into 1-D vector
  obs$cohort_peaks = as.numeric(peaks$peaks) # turn matrix into 1-D vector

  # ---- SIMULATED
  Lt   = .relative_age(date, birth, lifespan=lifespan)

  if(method=="rsa") {

    out = rsa(Linf=Linf, k=k, lifespan=lifespan, bins=bins, marks=marks,
              ASP=peaks$ASP, obs=obs, Lt=Lt, n=n, parallel=parallel, control=control)

  } else {

    out = NULL
    # out = calibrate(par, fn, gr, method=method, lower, upper, phases,
    #                 control, obs, Lt)

  }

  .t0 = Sys.time() - .t0
  out = c(out, list(elapsed=.t0))

  return(out)

}


ASP = function(x, return.peaks=FALSE) {
  peaks = apply(x, 2, .rlecount)
  peaks = peaks + max(peaks)*(col(peaks)-1)*(peaks!=0)
  mpeak = tapply(as.numeric(x), INDEX = as.numeric(peaks), FUN=max)
  if(isTRUE(return.peaks)) return(list(ASP=sum(mpeak), peaks=peaks))
  return(sum(mpeak))
}

# objective function

ESP = function(obs, sim, flag.method="first", ...) {

  if(nrow(sim)==0) return(0)
  # combine sim and obs
  sim = merge(sim, obs, sort=FALSE)
  esp = .ESP(x=sim, flag.method=flag.method)

  return(esp)
}

.ESP = function(x, flag.method) {
  if(!is.null(flag.method)) {
    ind = x$cohort_peaks==0
    neg = sum(x$Fs[ind], na.rm=TRUE)
    pos = sum(flag_out(x[!ind, ], flag.method=flag.method))
    esp = pos + neg
  } else {
    esp = sum(x$Fs, na.rm=TRUE)
  }
  return(esp)
}

rsa = function(Linf, k, lifespan, bins, marks, ASP, obs, Lt, n, parallel=FALSE, control) {

  # brute-force estimation squeleton
  if(length(Linf)==1) n[1] = 1
  Linfs = seq(from=Linf[1], to=Linf[2], length.out=n[1])
  ks    = seq(from=k[1], to=k[2], length.out=n[2])
  par = expand.grid(Linf=Linfs, k=ks)
  par$t0 = NA
  par$esp = NA

  pb = txtProgressBar(min=0, max=nrow(par), style=3, title = "RSA")
  opts = list(progress=function(n) setTxtProgressBar(pb, n))

  # setting para
  if(parallel) {

    if(isTRUE(control$init_cluster)) {
      cl = makeCluster(control$ncores, type="SOCK")
      registerDoSNOW(cl)
    }

  } else {

    registerDoSEQ()

  }

  .export = c("Lt", "lifespan", "bins", "marks", "obs")

  res = foreach(i=seq_len(nrow(par)), .combine=rbind, .options.snow=opts,
                .verbose=FALSE, .inorder=FALSE) %dopar% {

                  k    = par$k[i]
                  Linf = par$Linf[i]
                  t0   = .find_t0_dt(k=k, Linf=Linf, Lt=Lt, lifespan=lifespan, bins=bins, marks=marks, obs=obs)
                  if(!parallel) setTxtProgressBar(pb, i)
                  c(i, t0, attr(t0, "ESP"))

                }

  if(parallel & isTRUE(control$init_cluster)) stopCluster(cl)

  res = res[order(res[,1]), ]

  par$t0  = res[, 2]
  par$esp = res[, 3]

  opt = which.max(par$esp)
  opar = list(k=par$k[opt], Linf=par$Linf[opt])
  Rn = 0.1*10^(par$esp/ASP)
  rsa = list(Linf=Linfs, k=ks,
             Rn=matrix(Rn, nrow=length(Linfs), ncol=length(ks)))
  out = list(par=opar, Rn=max(Rn, na.rm=TRUE), rsa=rsa)

  return(out)

}

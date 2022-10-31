
#' @export
ELEFAN1 = function(x, ...) {
  UseMethod("ELEFAN1")
}

#' @export
ELEFAN1.len_freq = function(x, Linf, k, lifespan, birth=NULL, method="rsa",
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
  control = .elefan_control_defaults(control)

  if(is.null(birth)) {
    birth = "01-apr"
    if(is.null(control$fixed.birth)) control$fixed.birth = FALSE
  } else {
    if(is.null(control$fixed.birth)) control$fixed.birth = TRUE
  }

  .t0 = Sys.time()

  control$ncohorts = length(birth)

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
  relage   = .relative_age(date, birth, lifespan=lifespan, asp=peaks$asp_by_date)
  Lt = relage$Lt
  births = relage$birth

  if(method=="rsa") {

    out = rsa(Linf=Linf, k=k, lifespan=lifespan, bins=bins, marks=marks, births=births,
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
  ydate = col(peaks)
  mpeak = tapply(as.numeric(x), INDEX = as.numeric(peaks), FUN=max)
  ypeak = tapply(as.numeric(x), INDEX = as.numeric(ydate), FUN=max)
  if(isTRUE(return.peaks))
    return(list(ASP=sum(mpeak), peaks=peaks, asp_by_date=ypeak))
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

  if(flag.method=="noflag") flag.method = NULL

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

rsa = function(Linf, k, lifespan, bins, marks, births, ASP, obs, Lt, n, parallel=FALSE, control) {

  # brute-force estimation squeleton
  if(length(Linf)==1) n[1] = 1
  Linfs = seq(from=Linf[1], to=Linf[2], length.out=n[1])
  ks    = seq(from=k[1], to=k[2], length.out=n[2])
  if(any(k==0)) k[k==0] = 1e-3
  par = expand.grid(Linf=Linfs, k=ks)
  par$t0 = NA
  par$esp = NA

  ncoh = seq_len(diff(range(Lt$cohort, na.rm=TRUE)) + 1)

  ncohorts    = control$ncohorts
  flag.method = control$flag.method
  fixed.birth = control$fixed.birth

  pb = txtProgressBar(min=0, max=nrow(par), style=3, title = "RSA")
  opts = list(progress=function(n) setTxtProgressBar(pb, n))

  # setting parallel
  if(parallel) {

    if(isTRUE(control$init_cluster)) {
      cl = makeCluster(control$ncores, type="SOCK")
      registerDoSNOW(cl)
    }

  } else {

    registerDoSEQ()

  }

  .export = c("Lt", "lifespan", "ncohorts", "bins", "marks", "obs")

  res = foreach(i=seq_len(nrow(par)), .options.snow=opts,
                .verbose=FALSE, .inorder=FALSE) %dopar% {

      k    = par$k[i]
      Linf = par$Linf[i]
      t0   = .find_t0_dt(k=k, Linf=Linf, Lt=Lt, lifespan=lifespan,
                         bins=bins, marks=marks, obs=obs, ncohorts=ncohorts,
                         flag.method=flag.method, fixed.birth=fixed.birth)
      if(!parallel) setTxtProgressBar(pb, i)
      t0$i = i
      t0
      # c(i, t0, attr(t0, "ESP"), attr(t0, "range"), attr(t0, "bc_t0"), attr(t0, "bc_ESP"))

      }

  if(parallel & isTRUE(control$init_cluster)) stopCluster(cl)

  # sort results
  ind = order(sapply(res, FUN="[[", i="i"))
  res = res[ind]

  xres = sapply(res, FUN=function(x) c(x$t0, attr(x$t0, "ESP")))

  par$t0  = xres[1, ]
  par$esp = xres[2, ]

  opt = which.max(par$esp)
  opar = list(Linf=par$Linf[opt], k=par$k[opt], t0=par$t0[opt])
  Rn = 0.1*10^(par$esp/ASP)
  phi = outer(ks, Linfs, FUN=phi)

  ind_bc   = sapply(res, FUN=function(x) !is.null(x$bycohort))
  bc_esp   = do.call(c, lapply(res[ind_bc], FUN="[[", i="esp"))
  bycohort = do.call(rbind, lapply(res, FUN="[[", i="bycohort"))

  bycohort$ID = phi(bycohort$k, bycohort$Linf)
  colnames(bycohort)[4] = "phi"

  bcoh = .best_by_cohort(bycohort=bycohort, t0=opar$t0, births=births, Rn=0.1*10^(bc_esp/ASP))
  bcoh$ASP = NULL

  rsa = list(Linf = Linfs,
             k    = ks,
             t0   = matrix(par$t0, nrow=length(Linfs), ncol=length(ks)),
             phi  = phi,
             Rn   = matrix(Rn, nrow=length(Linfs), ncol=length(ks)),
             par  = par,
             bycohort = bycohort)

  out = list(par=opar, Rn=max(Rn, na.rm=TRUE), ASP=ASP, rsa=rsa, phi=phi[opt],
             bycohort=bcoh)

  return(out)

}

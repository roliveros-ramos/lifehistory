
#' @export
ELEFAN1 = function(x, ...) {
  UseMethod("ELEFAN1")
}

#' @export
ELEFAN1.len_freq = function(x, lifespan, Linf=NULL, k=NULL, sp=NULL, birth=NULL, method="arsa",
                            parallel=TRUE, control=list()) {

  control = .elefan_control_defaults(control)

  if(is.null(Linf) | is.null(k)) {
    if(method=="rsa") stop("Must provide 'Linf' and 'k' for method='rsa'")
    ng = .start_grid(x=x, MA=control$MA)
    if(is.null(Linf)) Linf = ng$Linf
    if(is.null(k))    k    = ng$k
  }

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

  method = match.arg(method, choices = c("rsa", "arsa"))
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

  Linf = .check_Linf(Linf, method=method)
  k    = .check_k(k, method=method)

  # ---- OBSERVED
  rcounts  = lf_restructure(x, MA=control$MA, method=control$method)
  peaks    = ASP(rcounts, return.peaks = TRUE)

  # new form of observed data
  obs = expand.grid(length=marks, date=date)
  obs$Fs = as.numeric(rcounts)  # turn matrix into 1-D vector
  maxFs = matrix(peaks$asp_by_date, nrow=nrow(rcounts), ncol=ncol(rcounts), byrow=TRUE)
  obs$maxFs = as.numeric(maxFs) # turn matrix into 1-D vector
  obs$cohort_peaks = as.numeric(peaks$peaks) # turn matrix into 1-D vector

  # ---- SIMULATED
  relage   = .relative_age(date, birth, lifespan=lifespan, asp=peaks$asp_by_date)
  Lt = relage$Lt
  births = relage$birth

  out = NULL

  if(method %in% c("rsa", "arsa")) {

    if(method=="rsa") {
      tmp = rsa(Linf=Linf, k=k, lifespan=lifespan, bins=bins, marks=marks, births=births,
                ASP=peaks$ASP, obs=obs, Lt=Lt, n=n, parallel=parallel, control=control)
    }

    if(method=="arsa") {
      tmp = arsa(Linf=Linf, k=k, lifespan=lifespan, bins=bins, marks=marks, births=births,
                ASP=peaks$ASP, obs=obs, Lt=Lt, n=n, parallel=parallel, control=control)
    }

    opt  = which.max(tmp$global$ESP) # par=global
    opar = list(Linf=tmp$global$Linf[opt], k=tmp$global$k[opt], t0=tmp$global$t0[opt])

    bcoh = .best_by_cohort(bycohort=tmp$bycohort, par=opar, births=births, Rn=tmp$bp_Rn)

    out = list(par=opar, lifespan=lifespan, births=births,
               Rn=max(tmp$Rn, na.rm=TRUE), ASP=peaks$ASP,
               phi=phi(k=opar$k, Linf=opar$Linf),
               rsa=tmp, bycohort=bcoh)

  } else {
    out = NULL
    # out = calibrate(par, fn, gr, method=method, lower, upper, phases,
    #                 control, obs, Lt)
  }

  # HERE to calculate best: global, bycohort, etc.

  .t0 = Sys.time() - .t0
  out = c(out, list(elapsed=.t0), list(control=control))

  return(out)

}


ASP = function(x, return.peaks=FALSE) {
  peaks = apply(x, 2, .rlecount)
  peaks = peaks + max(peaks)*(col(peaks)-1)*(peaks!=0)
  ydate = col(peaks)
  xsize = row(peaks)
  mpeak = tapply(as.numeric(x), INDEX = as.numeric(peaks), FUN=max)
  ypeak = apply(x, 2, max, na.rm=TRUE)
  npeaks = apply(peaks, 2, FUN=function(x) length(unique(x[x!=0])))

  xpeak = tapply(as.numeric(x), INDEX = as.numeric(peaks), FUN=which.max)
  tpeak = tapply(as.numeric(ydate), INDEX = as.numeric(peaks), FUN=min)
  speak = tapply(as.numeric(xsize), INDEX = as.numeric(peaks), FUN=identity)
  modes = sapply(seq_along(xpeak), FUN=function(i) speak[[i]][xpeak[i]])
  out = data.frame(cohort = as.numeric(names(tpeak)),
                   time   = tpeak,
                   row    = modes,
                   rmin   = sapply(speak, min),
                   rmax   = sapply(speak, max),
                   peak   = mpeak)
  if(isTRUE(return.peaks))
    return(list(ASP=sum(mpeak), peaks=peaks, asp_by_date=ypeak, npeaks=npeaks,
                data=out[-1, ]))
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

rsa = function(Linf, k, lifespan, bins, marks, births, ASP, obs, Lt, n, parallel=FALSE,
               use_grid=FALSE, control) {

  # brute-force estimation squeleton
  if(!isTRUE(use_grid)) {
    if(length(Linf)==1) {
      n[1] = 1
      Linf = c(Linf, Linf)
    }
    Linfs = seq(from=Linf[1], to=Linf[2], length.out=n[1])
    ks    = seq(from=k[1], to=k[2], length.out=n[2])
  } else {
    Linfs = Linf
    ks = k
  }
  if(any(ks==0)) ks[ks==0] = 1e-3
  par = expand.grid(Linf=Linfs, k=ks)
  par$t0 = NA
  par$ESP = NA

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
      on.exit(stopCluster(cl))
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
      }

  # sort results
  ind = order(sapply(res, FUN="[[", i="i"))
  res = res[ind]

  # add global t0 and ESP
  # xres = sapply(res, FUN=function(x) c(x$t0, attr(x$t0, "ESP")))
  # par$t0  = xres[1, ]
  # par$esp = xres[2, ]
  par$t0  = sapply(res, FUN=function(x) x$t0)
  par$phi = phi(k=par$k, Linf=par$Linf)
  par$ESP = sapply(res, FUN=function(x) attr(x$t0, "ESP"))
  par$Rn  = 0.1*10^(par$ESP/ASP)

  ind_bc   = sapply(res, FUN=function(x) !is.null(x$bycohort))
  bp_esp   = do.call(c, lapply(res[ind_bc], FUN="[[", i="esp"))
  bp_Rn    = 0.1*10^(bp_esp/ASP)
  bycohort = do.call(rbind, lapply(res, FUN="[[", i="bycohort"))

  bycohort$ID = phi(bycohort$k, bycohort$Linf)
  colnames(bycohort)[4] = "phi"

  rsa = list(Linf     = Linfs,
             k        = ks,
             t0       = matrix(par$t0, nrow=length(Linfs), ncol=length(ks)),
             Rn       = matrix(par$Rn, nrow=length(Linfs), ncol=length(ks)),
             phi      = matrix(par$phi, nrow=length(Linfs), ncol=length(ks)),
             bp_Rn    = bp_Rn, # more detail than Rn, more t0s.
             bycohort = bycohort, # ASP by cohort for each parameter
             global   = par)

  return(rsa)

}

# adaptative RSA

arsa = function(Linf, k, lifespan, bins, marks, births, ASP, obs, Lt, n, parallel=FALSE, control) {

  if(parallel) {

    if(isTRUE(control$init_cluster)) {
      cl = makeCluster(control$ncores, type="SOCK")
      registerDoSNOW(cl)
      on.exit(stopCluster(cl))
      control$init_cluster = FALSE
    }

  } else {

    registerDoSEQ()

  }

  out = list()

  for(i in seq_len(control$iter)) {

    msg = sprintf("\nIteration %s: k=[%0.2f, %0.2f], Linf=[%0.2f, %0.2f]",
                  i, min(k), max(k), min(Linf), max(Linf))
    message(msg)
    out[[i]] = rsa(Linf=Linf, k=k, lifespan=lifespan, bins=bins, marks=marks, births=births,
                   ASP=ASP, obs=obs, Lt=Lt, n=n, parallel=parallel, use_grid=TRUE, control=control)

    tmp = .update_grid(out[[i]], i, n)
    Linf = tmp$Linf
    k    = tmp$k

  }

  out = .merge_rsa(out)

  return(out)

}


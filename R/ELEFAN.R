
#' @export
ELEFAN1 = function(x, ...) {
  UseMethod("ELEFAN1")
}

#' @export
ELEFAN1.len_freq = function(x, Linf, k, lifespan, birth="01-apr", n = 30,
                            control=list()) {

  .t0 = Sys.time()

  marks = x$marks
  bins  = x$bins
  date  = x$date
  x     = x$x

  control = .elefan_control_defaults(control)

  n = .check_n(n)

  Linf = .check_Linf(Linf)
  k    = .check_k(k)

  Linfs = seq(from=Linf[1], to=Linf[2], length.out=n[1])
  ks    = seq(from=k[1], to=k[2], length.out=n[2])

  # ---- OBSERVED

  rcounts  = lf_restructure(x, MA=control$MA, method=control$method,
                            control=control)
  peaks    = ASP(rcounts, return.peaks = TRUE)

  # new form of observed data
  obs = expand.grid(length=marks, date=date)
  obs$Fs = as.numeric(rcounts)  # turn matrix into 1-D vector
  obs$cohort_peaks = as.numeric(peaks$peaks) # turn matrix into 1-D vector
  # obs = as.data.table(obs)
  # simulation squeleton
  par = expand.grid(Linf=Linfs, k=ks)
  par$t0 = NA
  par$esp = NA

  # ---- SIMULATED
  Lt   = .relative_age(date, birth, lifespan=lifespan)

  pb = txtProgressBar(min=0, max=nrow(par), style=3)

  for(i in seq_len(nrow(par))) {
    k    = par$k[i]
    Linf = par$Linf[i]
    t0   = .find_t0(k, Linf, Lt, lifespan, bins, obs)
    sim = cbind(Lt, length=.VB(Lt$age, Linf=Linf, k=k, t0=t0))
    # data.table experiment
    # sim = copy(Lt)
    # sim[, length := .VB(age, Linf=Linf, k=k, t0=t0)]
    par$esp[i] = ESP(obs, sim, bins, marks)
    par$t0[i] = t0
    setTxtProgressBar(pb, i)
  }

  .t0 = Sys.time() - .t0

  opt = which.max(par$esp)
  opar = list(k=par$k[opt], Linf=par$Linf[opt])
  Rn = 0.1*10^(par$esp/peaks$ASP)
  rsa = list(Linf=Linfs, k=ks,
             Rn=matrix(Rn, nrow=length(Linfs), ncol=length(ks)))
  out = list(par=opar, Rn=Rn, rsa=rsa, elapsed=.t0)
  class(out) = "elefan"
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

ESP = function(obs, sim, bins, marks, flag.method="first", ...) {

  # transform length to length-class
  sim$length = marks[cut(sim$length, breaks=bins, right=FALSE, labels=FALSE)]
  # remove length outside sampling
  sim = sim[!is.na(sim$length), ]
  if(nrow(sim)==0) return(0)
  # combine sim and obs
  sim = merge(sim, obs, sort=FALSE)

  if(!is.null(flag.method)) {
    ind = sim$cohort_peaks==0
    neg = sum(sim$Fs[ind], na.rm=TRUE)
    pos = sum(flag_out(sim[!ind, ], flag.method=flag.method))
    esp = pos + neg
  } else {
    esp = sum(sim$Fs, na.rm=TRUE)
  }

  return(esp)
}





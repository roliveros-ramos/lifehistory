
# Counting peaks detected after restructuring -----------------------------

#' Counting peaks detected after restructuring
#' @param x A 'len_freq' object created with lf_create.
#' @param N The maximum MA tested is 2N-1.
#'
#' @export
lf_modes = function(x, N=7) {
  out = sapply(1:N, FUN=.npeaks, x=x$x)
  colnames(out) = sprintf("MA=%d", 2*(1:N)+1)
  class(out) = c("len_freq_modes", class(out))
  return(out)
}

lf_analysis = function(x, MA) {
  MA = as.integer(MA)
  if(MA%%2 == 0) stop("MA must be an odd integer.")
  rx = lf_restructure(x$x, MA = MA)
  asp = ASP(rx, return.peaks = TRUE)
  dat = asp$data
  dat$mode = x$marks[dat$row]
  dat$date = x$date[dat$time]
  dat$size_min = x$marks[dat$rmin]
  dat$size_max = x$marks[dat$rmax]
  dat$max = unlist(tapply(dat$mode, INDEX = dat$date, FUN=function(x) x==max(x)))
  dat$max = as.numeric(dat$max)
  dat = dat[, c("date", "mode", "size_min", "size_max", "peak", "max")]

  lmax = tapply(dat$mode[dat$max==1], INDEX = year(dat$date)[dat$max==1], FUN = max)

  attr(dat, "lmax") = lmax
  return(dat)
}

#' @export
plot.len_freq_modes = function(x, ...) {
  boxplot(x, ylim=c(0, 1.1*max(x)), las=1)
  return(invisible(TRUE))
}

.npeaks = function(n, x) {
  rx = lf_restructure(x, MA = 2*n+1)
  asp = ASP(rx, return.peaks = TRUE)$npeaks
  return(asp)
}


# Restructure length frequencies ------------------------------------------


lf_restructure = function(x, MA=5, method="GP1997", control=list(), ...) {

  # add check/warning for Lmax

  MA = as.integer(MA)
  if(MA%%2 == 0) stop("MA must be an odd integer.")

  .lf = switch(method,
               GP1997 = .lf_restructure_GP1997,
               BP1986 = .lf_restructure_BP1986,
               stop("Undefined method.")
  )

  # default options
  lowfreq = if(is.null(control$lowfreq)) FALSE else control$lowfreq
  trim    = if(is.null(control$trim))    TRUE else control$trim
  before  = if(is.null(control$before))  FALSE else control$before

  out = apply(x, MARGIN=2, FUN=.lf, MA=MA, lowfreq=lowfreq, trim=trim, before=before)
  return(out)

}


# Brey & Pauly 1986
.lf_restructure_BP1986 = function(x, MA, lowfreq=TRUE, trim=FALSE, before=TRUE) {

  after = !before
  ind = .trimZeros(x)
  Fi = if(isTRUE(trim)) x[ind] else x

  m = (MA-1)/2
  n = length(Fi)
  xFi = c(rep(0, m), Fi, rep(0, m)) # add zero padding
  # ----- STEP 1: moving average (MAi)
  MAi = filter(xFi, filter=rep(1, MA)/MA)
  MAi = tail(head(MAi, -m), -m)
  # ----- STEP 2: transformed length-frequency (fi)
  AFi = Fi/MAi
  # ----- STEP 3: rescaled lf
  fi = AFi/mean(AFi, na.rm=TRUE) - 1
  # ----- STEP 4: rescaling positives
  # number of zeros (nz)
  nz  = filter(xFi==0, rep(1, MA))
  nz  = tail(head(nz, -m), -m)
  nz  = ifelse(fi>0, nz, 0) # only corrects positives
  fi = fi*(1/2^nz)
  # optional correction for frequency
  if(isTRUE(lowfreq)) {
    fc = ifelse(fi>0, sqrt(1 + 2/Fi^2), 1)
    fi = fi/fc
  }
  # ----- STEP 5: correction of negatives (sum AFTER setting zeros)
  if(isTRUE(before)) fi[which(Fi==0)] = 0 # -1 to 0 for Fi==0
  sump = sum(fi*(fi>0))
  sumn = -sum(fi*(fi<0))
  if(isTRUE(after)) fi[which(Fi==0)] = 0
  CF = ifelse(fi<0, sump/sumn, 1)
  fi = fi*CF
  # ----- STEP 6: correction of tails

  if(isTRUE(trim)) {
    out = numeric(length(x))
    out[ind] = fi
    fi = out
  }

  if(fi[n]<0) fi[n] = 0
  if(fi[n-1]<0) fi[n-1] = 0.5*fi[n-1]

  return(fi)
}

# Gayanilo & Pauly 1997
.lf_restructure_GP1997 = function(x, MA, lowfreq=FALSE, trim=TRUE, before=FALSE) {

  after = !before
  ind = .trimZeros(x)
  Fi = if(isTRUE(trim)) x[ind] else x
  m = (MA-1)/2
  n = length(Fi)
  xFi = c(rep(0, m), Fi, rep(0, m)) # add zero padding
  # ----- STEP A: moving average (MAi)
  MAi = filter(xFi, filter=rep(1, MA)/MA)
  MAi = tail(head(MAi, -m), -m)
  # ----- STEP B: transformed length-frequency (fi)
  AFi = Fi/MAi
  AFi[is.na(AFi) & Fi==0] = 0
  mp  = mean(AFi, na.rm=TRUE)
  # ----- STEP C: rescaled lf
  fi = AFi/mp - 1
  # ----- STEP D: Isolated peaks, count zeros (nz)
  nz  = filter(xFi==0, rep(1, MA))
  nz  = tail(head(nz, -m), -m)
  nz  = ifelse(fi>0, nz, 0) # only corrects positives
  # ----- STEP E: De-emphasizing positive
  fi = fi*(0.5^nz)
  # optional correction for frequency
  if(isTRUE(lowfreq)) {
    fc = ifelse(fi>0, sqrt(1 + 2/Fi^2), 1)
    fi = fi/fc
  }

  # correction of Lmax: this is a global correction, but for the date?
  n = length(fi)
  if(fi[n]<0) fi[n] = 0
  if(fi[n-1]<0) fi[n-1] = 0.5*fi[n-1]

  # ----- STEP F: correction of negatives (sum BEFORE setting 0s)
  if(isTRUE(before)) fi[which(fi==-1)] = 0
  sump = sum(fi*(fi>0))
  sumn = -sum(fi*(fi<0))
  if(isTRUE(after)) fi[which(fi==-1)] = 0 # -1 to 0 for Fi==0
  CF = ifelse(fi<0, sump/sumn, 1)
  fi = fi*CF

  if(isTRUE(trim)) {
    out = numeric(length(x))
    out[ind] = fi
    fi = out
  }

  return(fi)
}


# ELEFAN internals --------------------------------------------------------

.check_n = function(n) {
  if(length(n)>2) stop("A maximum of two 'n' values (n_Linf, n_k) can be provided.")
  if(length(n)==1) n = c(n, n)
  if(any(n<0)) stop("All 'n' must be positive integers.")
  n = as.integer(n)
  n = n + 1 - n%%2 # always odd
  return(n)
}


.check_Linf = function(Linf, method=NULL) {

  if(method=="arsa") {
    if(any(Linf<=0)) stop("Only positive values for Linf are valid.")
    return(sort(Linf))
  }

  if(missing(Linf)) stop("You must provide at least one 'Linf' value.")
  if(length(Linf)==0) stop("You must provide at least one 'Linf' value.")
  if(length(Linf)==1) {
    message("Using fix Linf, profiling over 'k' curvature parameter.")
  }
  if(length(Linf)>2) {
    stop("A range (min, max) for 'Linf' search must be provided.")
  }
  if(any(Linf<=0)) stop("Only positive values for Linf are valid.")
  if(Linf[1]>Linf[2]) stop("Increasing values for Linf are expected.")
  if(Linf[1]==Linf[2]) {
    message("Identical minimum and maximum search values for Linf have been provided.")
    message("Using fix Linf, profiling over 'k' curvature parameter.")
    Linf = Linf[1]
  }

  return(Linf)

}

.check_k = function(k, method=NULL) {

  if(method=="arsa") {
    if(any(k<0)) stop("Only positive values for k are valid.")
    return(sort(k))
  }

  if(missing(k)) stop("You must provide at least two 'k' value.")
  if(length(k)<2) stop("You must provide at least two 'k' value.")
  if(length(k)>2) {
    stop("A range (min, max) for 'k' search must be provided.")
  }
  if(any(k<0)) stop("Only positive values for k are valid.")
  if(k[1]>k[2]) stop("Increasing values for k are expected.")
  if(k[1]==k[2]) stop("You must provide at least two different 'k' values.")

  return(k)

}

.elefan_control_defaults = function(control) {
  if(is.null(control$MA)) control$MA = 5
  if(is.null(control$flag.method)) control$flag.method = "first"
  if(is.null(control$method)) control$method = "GP1997"
  # default options for length restructure
  if(is.null(control$lowfreq)) control$lowfreq = FALSE
  if(is.null(control$trim))    control$trim = FALSE
  if(is.null(control$before))  control$before = FALSE
  if(is.null(control$ncores))  control$ncores = detectCores(logical = FALSE)
  if(is.na(control$ncores))    control$ncores = 1
  if(is.null(control$iter))    control$iter = 3
  if(is.null(control$n))       control$n = 10
  control$n = .check_n(control$n)
  if(is.null(control$init_cluster)) control$init_cluster=TRUE
  if(is.null(control$fixed.birth)) control$fixed.birth = FALSE
  return(control)
}


# ELEFAN methods ----------------------------------------------------------

#' @export
plot.elefan = function(x, what="Rn", zlim=NULL, p=1.42, col=NULL, ...) {

  x = x$rsa
  y = x[[what]]

  thr = 0.1
  zz = if(what=="Rn") c(0,1) else NULL
  if(is.null(zlim)) zlim=range(c(y, zz), na.rm=TRUE)
  if(what=="Rn" & is.null(col)) col = elefanPalette(n=256, p=p)
  if(what=="t0" & is.null(col)) col = colorful::divergencePalette(zlim=zlim, center=0, p=p, symmetric = TRUE)
  image.plot(x$Linf, x$k, y, col=col, zlim=zlim,
             xlab = "Asymptotic length (Linf, cm)",
             ylab = "Curvature parameter (k, year-1)",
             main = sprintf("Response Surface Analysis (%s)", what))
  return(invisible(NULL))
}

#' @export
print.elefan = function(x) {
  cat("Best growth parameters by cohort:\n")
  print(x$bycohort)
  cat("\nOverall best growth parameters:\n")
  print(unlist(x$par))
  cat(sprintf("\nAssumed lifespan: %d years.\n", x$lifespan))
}


elefan.col = c("white", "cornsilk", "cornsilk",
               "khaki1", "orange", "darkorange",
               "red", "red3", "firebrick4",
               "violetred4")

elefanPalette = function(n, p) colorRampPalette(elefan.col, bias=p)(n)




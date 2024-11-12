#' Get priors for lifehistory parameters
#'
#' @param taxon The taxon you want to the priors.
#' @param par The parameters you want, possible options are "Linf", "k", "Winf",
#' "tmax", "tm", "M", "Lm" and "temperature".
#'
#' @return A list containing the mean and covariance for the given parameters.
#' @export
get_prior = function(taxon, par=NULL) {

  nm = c("Linf", "k", "Winf", "tmax", "tm", "M", "Lm", "temperature")
  if(is.null(par)) par = nm
  par = match.arg(par, choices=nm, several.ok = TRUE)

  x = get_classification(taxon, db = "fishlife")
  x = x[c("class", "order", "family", "genus", "species")]
  if(!is.na(x$species)) x$species = tail(strsplit(x$species, split=" ")[[1]], 1)
  x = unlist(x)
  x[is.na(x)] = "predictive"
  .ancestors = function(n, x) {
    x = rev(x)
    x[seq_len(n)] = "predictive"
    out = paste(rev(x), collapse="_")
    return(out)
  }
  out = lapply(5:0, FUN=.ancestors, x=x)
  names(out) = c("default", names(x))

  .get_mean = function(x) {
    lmu = rep(NA, length(nm))
    ind = which(rownames(priors$mean)==x)
    if(length(ind)==1) lmu = priors$mean[x, ]
    names(lmu) = nm
    lmu = lmu[par]
    return(lmu)
  }

  .get_sigma = function(x) {
    lsigma = matrix(NA, nrow=length(nm), ncol=length(nm))
    ind = which(rownames(priors$mean)==x)
    if(length(ind)==1) lsigma = priors$cov[x, ,]
    rownames(lsigma) = nm
    colnames(lsigma) = nm
    lsigma = lsigma[par, par]
    return(lsigma)
  }

  lmu = lapply(out, FUN=.get_mean)
  lsigma = lapply(out, FUN=.get_sigma)

  x = mapply(FUN=list, lmu, lsigma, out, SIMPLIFY = FALSE)
  x = lapply(x, FUN=setNames, nm=c("mean", "covariance", "taxon"))
  return(x)
}


# Auxiliar ----------------------------------------------------------------

.new_par = function(x, y, n, r=1) {
  qk = cumsum(x)/sum(x)
  fk = splinefun(x=qk, y=y)
  fpk = splinefun(x=qk, y=x)
  nq = seq(from=min(qk), to=max(qk), length.out=n)
  nk = fk(nq)
  pk = fpk(nq)
  nk = unique(round(nk/r, 0)*r)
  nk = base::setdiff(nk, y)
  return(nk)
}


.start_grid = function(x, MA) {

  dat = lf_analysis(x, MA=MA)
  lmax = attr(dat, "lmax")
  lmax = round(quantile(lmax, prob=c(0.05, 0.95))*c(0.9, 1.1))

  Linf = log(seq(lmax[1], lmax[2], by=0.01))
  k    = log(seq(0, 4, by=0.005))
  par = expand.grid(Linf=Linf, k=k)
  par$prob = dmvnorm(x=par[, c("Linf", "k")],
                     mean=default_prior$mean,
                     sigma=default_prior$covariance)
  prob = matrix(par$prob, nrow=length(Linf), ncol=length(k))
  prob = prob/sum(prob)

  pLinf = rowSums(prob)
  pk    = colSums(prob)

  nLinf = .new_par(x=pLinf, y=exp(Linf), n=30, r=0.5)
  nk    = .new_par(x=pk, y=exp(k), n=30, r=0.1)

  return(list(Linf=nLinf, k=nk))

}



.update_grid = function(x, i, n) {

  steps = list(k = c(0.1, 0.05, 0.01, 0.001),
               Linf = c(0.5, 0.1, 0.05, 0.01))

  par = x$global
  obj = x$Rn
  Linf = x$Linf
  k = x$k

  pLinf = rowMeans(obj)
  pk    = colMeans(obj)

  nLinf = .new_par(x=pLinf, y=Linf, n=n[1], r=steps$Linf[i+1])
  nk    = .new_par(x=pk, y=k, n=n[2], r=steps$k[i+1])

  return(list(Linf=nLinf, k=nk))

}

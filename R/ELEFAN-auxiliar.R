

# Auxiliar ----------------------------------------------------------------

.relative_age = function(date, birth, lifespan) {
  dbirth = yday(parse_date_time(birth, order=c("dm","md")))
  # t_anchor = dbirth/365
  start = min(date) - years(lifespan)
  yday(start) = min(dbirth)
  end = max(date)
  years = seq(from=year(start), to=year(end))
  nd = rep(rep(2000, length=length(years)), each=length(birth))
  # bs is the assumed birthdate of all cohorts (t_anchor)
  bd = date(sprintf("%s-01-01", nd))
  yday(bd) = dbirth
  year(bd) = rep(years, each=length(birth))
  age = outer(date, bd, difftime, units="days")/365.25
  Lt = data.frame(date=date, cohort=as.numeric(col(age)), age=as.numeric(age))
  Lt = Lt[Lt$age<=lifespan & Lt$age>0, ]
  # Lt = as.data.table(Lt)
  return(Lt)
}

.trimZeros = function(x) {
  rr = range(which(x!=0))
  ind = seq(from=rr[1], to=rr[2])
  return(ind)
}

.rlecount = function(x) {
  r = rle(ifelse(x>0, 1, 0))
  r$values[r$values==1] = seq_len(sum(r$values))
  return(inverse.rle(r))
}


.calculate_bins = function(marks) {
  bw = diff(marks) # bin width (should allow for uneven bin sizes)
  b_lower = marks - (c(bw[1], bw)/2) # upper bin limit
  b_upper = marks + (c(bw, bw[length(bw)])/2) # lower bin limit
  bins = sort(unique(c(b_lower, b_upper)))
  return(bins)
}


flag_out = function(x, flag.method="first") {
  if(nrow(x)==0) return(0)
  return(sapply(split(x, f=x$cohort_peaks), .flag, flag.method=flag.method))
}

.flag = function(x, flag.method) {
  out = switch(flag.method,
               "first" = x$Fs[which.min(x$cohort)],
               "last"  = x$Fs[which.max(x$cohort)],
               "max"   = max(x$Fs, na.rm=TRUE),
               "mean"  = mean(x$Fs, na.rm=TRUE))
  return(out)
}

.flag2 = function(x, flag.method) {
  out = switch(flag.method,
               "first" = x[which.min(cohort), .(Fs)],
               "last"  = x[which.max(cohort), .(Fs)],
               "max"   = x[, max(Fs, na.rm=TRUE)],
               "mean"  = x[, mean(Fs, na.rm=TRUE)])
  return(out)
}



.find_t0 = function(k, Linf, Lt, lifespan, bins, marks, obs, by=NULL, range=FALSE) {

  if(is.null(by)) {
    rr = .find_t0(k=k, Linf=Linf, Lt=Lt, lifespan=lifespan, bins=bins,
                  marks=marks, obs=obs, by=0.05, range=TRUE)
    by = 0.001
  } else {
    rr = c(-lifespan, 0)
  }

  t0s = seq(from=rr[1], to=rr[2], by=by)
  esp = numeric(length(t0s))

  for(i in seq_along(t0s)) {
    sim = cbind(Lt, length=.VB(Lt$age, Linf=Linf, k=k, t0=t0s[i])) #data.table
    esp[i] = ESP(obs, sim, bins, marks)
  }
  t0 = t0s[which.max(esp)]

  if(isTRUE(range)) {
    esp = esp/max(esp, na.rm=TRUE)
    return(range(t0s[esp>0.8]))
  }

  attr(t0, "ESP") = max(esp)
  return(t0)
}


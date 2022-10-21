

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

# .flag2 = function(x, flag.method) {
#   out = switch(flag.method,
#                "first" = x[which.min(cohort), .(Fs)],
#                "last"  = x[which.max(cohort), .(Fs)],
#                "max"   = x[, max(Fs, na.rm=TRUE)],
#                "mean"  = x[, mean(Fs, na.rm=TRUE)])
#   return(out)
# }



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
    # transform length to length-class
    sim$length = marks[cut(sim$length, breaks=bins, right=FALSE, labels=FALSE)]
    # remove length outside sampling
    sim = sim[!is.na(sim$length), ]
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


.find_t0_dt = function(k, Linf, Lt, lifespan, bins, marks, obs, by=NULL, range=FALSE) {

  if(is.null(by)) {
    rr = .find_t0_dt(k=k, Linf=Linf, Lt=Lt, lifespan=lifespan, bins=bins,
                  marks=marks, obs=obs, by=1/100, range=TRUE)
    by = 1/365
  } else {
    rr = c(-1, 0)
  }

  # if(!is.null(attr(rr, "ESP"))) return(rr)

  t0s = seq(from=rr[1], to=rr[2], by=by)
  par = expand.grid(Linf=Linf, k=k, t0=t0s)
  par$ID = seq_len(nrow(par))
  sim = data.frame(ID = par$ID,
                   t0 = par$t0,
                   date=rep(Lt$date, nrow(par)),
                   cohort=rep(Lt$cohort, nrow(par)),
                   age=rep(Lt$age, nrow(par)))
  simL = .VB(sim$age, Linf, k, sim$t0)
  sim$length = marks[cut(simL, breaks=bins, labels=FALSE)]
  sim = sim[!is.na(sim$length), ]

  if(nrow(sim)==0) {

    little = sum(simL < bins[1], na.rm=TRUE)
    big    = sum(simL > tail(bins, 1), na.rm=TRUE)

    t0 = ifelse(little <= big, -1, 0)
    attr(t0, "ESP") = -99

  } else {
    sim = merge(sim, obs, sort=FALSE)

    xesp = sapply(split(sim, f=sim$ID, drop = FALSE), FUN=.ESP, flag.method="first")

    esp = rep(NA_real_, length(t0s))
    esp[as.numeric(names(xesp))] = xesp

    t0 = t0s[which.max(esp)]
    attr(t0, "ESP") = max(esp, na.rm=TRUE)
  }

  if(isTRUE(range)) {
    xrr = pmin(pmax(t0 + 0.5*c(-1,1)/12, -1), 0)
    return(xrr)
  }

  return(t0)

}


.find_t0_dt2 = function(k, Linf, Lt, lifespan, bins, marks, obs) {

  t0s = seq(from=-1, to=0, length.out=365)
  par = expand.grid(Linf=Linf, k=k, t0=t0s)
  par$ID = seq_len(nrow(par))
  sim = data.frame(ID = par$ID,
                   t0 = par$t0,
                   date=rep(Lt$date, nrow(par)),
                   cohort=rep(Lt$cohort, nrow(par)),
                   age=rep(Lt$age, nrow(par)))
  sim$length = marks[cut(.VB(sim$age, Linf, k, sim$t0), breaks=bins, labels=FALSE)]
  sim = sim[!is.na(sim$length), ]
  sim = merge(sim, obs, sort=FALSE)

  xesp = sapply(split(sim, f=sim$ID, drop = FALSE), FUN=.ESP, flag.method="first")

  esp = rep(NA_real_, length(t0s))
  esp[as.numeric(names(xesp))] = xesp

  t0 = t0s[which.max(esp)]
  attr(t0, "ESP") = max(esp, na.rm=TRUE)

  return(t0)
}

# new optimized version ---------------------------------------------------


ESP2 = function(obs, sim, bins, marks, flag.method="first", ...) {

  # transform length to length-class
  sim$length = marks[cut(sim$length, breaks=bins, right=FALSE, labels=FALSE)]
  # remove length outside sampling
  sim = sim[!is.na(sim$length), ]
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


# if(!is.null(flag.method)) {
#   ind = sim$cohort_peaks==0
#   neg = sum(sim$Fs[ind], na.rm=TRUE)
#   pos = sum(flag_out(sim[!ind, ], flag.method=flag.method))
#   esp = pos + neg
# } else {
#   esp = sum(sim$Fs, na.rm=TRUE)
# }


# replace marks with integers in obs (no marks needed)
# sim[, length := cut(length, breaks=bins, right=FALSE, labels=FALSE)]

ESP2 = function(obs, sim, bins, marks, flag.method="first", ...) {

  # transform length to length-class
  sim[, length := cut(length, breaks=bins, right=FALSE, labels=FALSE)]
  # remove length outside sampling
  sim = sim[!is.na(length)]
  sim = obs[sim, on=c("date", "length")]

  if(!is.null(flag.method)) {
    # ind = sim$cohort_peaks==0
    neg = sim[cohort_peaks==0, sum(Fs, na.rm=TRUE)]
    pos = sim[cohort_peaks!=0, .flag2(.SD, flag.method=flag.method), by=cohort_peaks]
    esp = pos[, sum(Fs)] + neg
  } else {
    esp = sim[, sum(Fs, na.rm=TRUE)]
  }

  return(esp)
}

.flag_dt = function(x, flag.method="first") {

  if(nrow(x)==0) return(0)
  if(!is.null(flag.method)) {
    print(class(x))

    neg = x[cohort_peaks==0, sum(Fs, na.rm=TRUE)]
    pos = x[cohort_peaks!=0, .flag2(.SD, flag.method=flag.method), by=cohort_peaks]
    esp = pos[, sum(Fs)] + neg
  } else {
    esp = x[, sum(Fs, na.rm=TRUE)]
  }

  return(esp)

}

.xflag = function(x, flag.method="first") {

  if(nrow(x)==0) return(0)

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





.flag2 = function(x, flag.method) {
  out = switch(flag.method,
               "first" = as.numeric(x[which.min(cohort), .(Fs)]),
               "last"  = as.numeric(x[which.max(cohort), .(Fs)]),
               "max"   = x[, max(Fs, na.rm=TRUE)],
               "mean"  = x[, mean(Fs, na.rm=TRUE)])
  return(out)
}

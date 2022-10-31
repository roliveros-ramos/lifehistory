

# Auxiliar ----------------------------------------------------------------

.relative_age = function(date, birth, lifespan, asp) {
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
  bd = bd[bd < end]
  births = data.frame(cohort=seq_along(bd), bd)
  age = outer(date, bd, difftime, units="days")/365.25
  Lt = data.frame(date=date, cohort=as.numeric(col(age)), age=as.numeric(age))
  Lt = Lt[Lt$age<=lifespan & Lt$age>0, ]
  ncounts = table(Lt$cohort)
  ncounts = data.frame(cohort=as.numeric(names(ncounts)), n=as.numeric(ncounts))
  births = merge(births, ncounts)

  asp_date = data.frame(date=date, asp=asp)
  asp_date = merge(Lt, asp_date)
  asp_date = aggregate(asp_date$asp, by=list(cohort=asp_date$cohort), FUN=sum)
  births = merge(births, asp_date)
  births$ASP = births$x
  births$x = NULL

  # Lt = as.data.table(Lt)
  return(list(Lt=Lt, birth=births))
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



.find_t0 = function(k, Linf, Lt, lifespan, bins, marks, obs, ncohort, by=NULL, range=FALSE) {

  if(is.null(by)) {
    rr = .find_t0(k=k, Linf=Linf, Lt=Lt, lifespan=lifespan, bins=bins,
                  marks=marks, obs=obs, by=0.05, range=TRUE)
    by = 0.001
  } else {
    rr = c(-0.5, 0.5)
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


.find_t0_dt = function(k, Linf, Lt, lifespan, bins, marks, obs, ncohorts,
                       flag.method, fixed.birth=FALSE, by=NULL, range=FALSE) {


  null_value = -99

  if(isTRUE(fixed.birth)) {

    rr = c(-1, 1)/12 #
    by = 1/365
    range = FALSE

  } else {

    if(is.null(by)) {
      rr = .find_t0_dt(k=k, Linf=Linf, Lt=Lt, lifespan=lifespan, bins=bins,
                       marks=marks, obs=obs, ncohorts=ncohorts, flag.method=flag.method,
                       fixed.birth=FALSE, by=1/100, range=TRUE)
      by = 1/365
    } else {
      rr = c(-0.5, 0.5)/ncohorts
    }

  }

  # if(!is.null(attr(rr, "ESP"))) return(rr)

  t0s = seq(from=rr[1], to=rr[2], by=by)
  par = expand.grid(Linf=Linf, k=k, t0=t0s)

  IDs = seq_len(nrow(par))
  coh = seq(from=min(Lt$cohort), to=max(Lt$cohort))
  bycohort = matrix(NA_real_, nrow=length(IDs), ncol=length(coh))

  par$ID = IDs
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

    t0 = ifelse(little <= big, -0.5, 0.5)/ncohorts
    attr(t0, "ESP") = null_value

    out = list(t0=t0, bycohort=NULL, esp=null_value)
    # attr(t0, "bc_ESP") = rep(null_value, length(coh))
    # attr(t0, "bc_t0") = rep(t0, length(coh))
    # attr(t0, "range") = rr

  } else {

    sim = merge(sim, obs, sort=FALSE)

    xesp = sapply(split(sim, f=sim$ID, drop = FALSE), FUN=.ESP, flag.method=flag.method)

    esp = rep(NA_real_, length(t0s))
    esp[as.numeric(names(xesp))] = xesp

    t0 = t0s[which.max(esp)]
    attr(t0, "ESP") = max(esp, na.rm=TRUE)

    if(!isTRUE(range)) {

      cx = tapply(sim$Fs, INDEX=sim[, c("ID", "cohort")], FUN=sum)
      xrow = match(rownames(cx), IDs)
      xcol = match(colnames(cx), coh)
      bycohort[xrow, xcol] = cx
      bycohort[is.na(bycohort)] = null_value
      bycohort = as.data.frame(bycohort)
      colnames(bycohort) = sprintf("cohort_%d", coh)
      bycohort = cbind(par, bycohort)

      # here, filter IDs within a range of t0 (global)
      # bc_esp = apply(bycohort, 2, max)
      # bc_t0s = t0s[apply(bycohort, 2, which.max)]

      # attr(t0, "bc_ESP") = bc_esp
      # attr(t0, "bc_t0") = bc_t0s
      # attr(t0, "range") = rr

      out = list(t0=t0, bycohort=bycohort, esp=esp)

    }

  }

  if(isTRUE(range)) {
    xrr = pmin(pmax(t0 + 0.5*c(-1,1)/12, -0.5/ncohorts), 0.5/ncohorts)
    return(xrr)
  }

  return(out)

}

.find_t0_dt_vector = function(k, Linf, Lt, lifespan, bins, marks, obs, ncohorts,
                       flag.method, fixed.birth=FALSE, by=NULL, range=FALSE) {


  null_value = -99

  if(isTRUE(fixed.birth)) {

    rr = c(-1, 1)/12 #
    by = 1/365
    range = FALSE

  } else {

    if(is.null(by)) {
      rr = .find_t0_dt_vector(k=k, Linf=Linf, Lt=Lt, lifespan=lifespan, bins=bins,
                       marks=marks, obs=obs, ncohorts=ncohorts, flag.method=flag.method,
                       fixed.birth=FALSE, by=1/100, range=TRUE)
      by = 1/365
    } else {
      rr = c(-0.5, 0.5)/ncohorts
    }

  }

  # if(!is.null(attr(rr, "ESP"))) return(rr)

  t0s = seq(from=rr[1], to=rr[2], by=by)
  par = expand.grid(Linf=Linf, k=k, t0=t0s)

  IDs = seq_len(nrow(par))
  coh = seq(from=min(Lt$cohort), to=max(Lt$cohort))
  bycohort = matrix(NA_real_, nrow=length(IDs), ncol=length(coh))

  par$ID = IDs
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

    t0 = ifelse(little <= big, -0.5, 0.5)/ncohorts
    attr(t0, "ESP") = null_value
    attr(t0, "bc_ESP") = rep(null_value, length(coh))
    attr(t0, "bc_t0") = rep(t0, length(coh))
    attr(t0, "range") = rr

  } else {

    sim = merge(sim, obs, sort=FALSE)

    xesp = sapply(split(sim, f=sim$ID, drop = FALSE), FUN=.ESP, flag.method=flag.method)

    esp = rep(NA_real_, length(t0s))
    esp[as.numeric(names(xesp))] = xesp

    t0 = t0s[which.max(esp)]
    attr(t0, "ESP") = max(esp, na.rm=TRUE)

    if(!isTRUE(range)) {

      cx = tapply(sim$Fs, INDEX=sim[, c("ID", "cohort")], FUN=sum)
      xrow = match(rownames(cx), IDs)
      xcol = match(colnames(cx), coh)
      bycohort[xrow, xcol] = cx
      bycohort[is.na(bycohort)] = null_value
      # here, filter IDs within a range of t0 (global)
      bc_esp = apply(bycohort, 2, max)
      bc_t0s = t0s[apply(bycohort, 2, which.max)]

      attr(t0, "bc_ESP") = bc_esp
      attr(t0, "bc_t0") = bc_t0s
      attr(t0, "range") = rr

    }

  }

  if(isTRUE(range)) {
    xrr = pmin(pmax(t0 + 0.5*c(-1,1)/12, -0.5/ncohorts), 0.5/ncohorts)
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

.best_by_cohort = function(bycohort, t0, births, Rn, p=0.5) {
  if(length(Rn)!=nrow(bycohort)) stop("global and by-cohort ESP dimension do not match.")
  if(p<0 | p>1) stop("Weighting factor 'p' must be between 0 and 1.")
  Rn = as.numeric(Rn)
  t1 = t0 + c(-1, 1)/12 # one month of variability for recruitment peak
  # bcoh = bycohort[bycohort$t0>=t1[1] & bycohort$t0<=t1[2], ]
  bcoh = bycohort
  xasp = matrix(births$ASP, nrow=nrow(bcoh), ncol=nrow(births), byrow=TRUE)
  bc_esp = 0.1*10^(as.matrix(bcoh[, -c(1:4)])/xasp)
  bc_esp = p*bc_esp + (1-p)*Rn
  bcoh = bcoh[apply(bc_esp, 2, which.max), ]
  cohs = gsub(colnames(bcoh)[-c(1:4)], pattern="cohort_", replacement="")
  bcoh = cbind(cohort=as.numeric(cohs), bcoh[, 1:4], ESP=diag(as.matrix(bcoh[, -c(1:4)])))
  rownames(bcoh) = NULL
  bcoh = merge(bcoh, births)
  bcoh$cohort = bcoh$bd
  bcoh$bd = NULL
  bcoh$ESP = NULL
  return(bcoh)
}



# Growth models -----------------------------------------------------------

phi = function(k, Linf) log10(k) + 2*log10(Linf)

# von Bertallanfy
.VB = function(x, Linf, k, t0) Linf*(1 - exp(-k*(x-t0)))

.VB_inv = function(x, Linf, k, t0) -(log(1-x/Linf) - k*t0)/k

#Quince-Boukal-Dieckmann

.QBD_immature = function(age, L0, a, b, c, beta, r=0, am=NULL, age0=0) {
  age = age - age0
  L = (L0^(b*(1-beta)) + c*(1-beta)/(a^(1-beta))*age)^(1/(b*(1-beta)))
  return(L)
}

.QBD_mature = function(age, L0, a, b, c, beta, r, am, age0=0) {
  age = age - age0
  am  = am -  age0
  R = 1/(1 + (1-beta)*r)
  H = c*(1-beta)/(a^(1-beta))
  L = ((R^(age-am))*(L0^(b*(1-beta)) + H*am) + (1 - R^(age-am))*R*H/(1-R))^(1/(b*(1-beta)))
  return(L)
}

.QBD = function(am, age, L0, a, b, c, beta, r, age0=0) {
  Lj = .QBD_immature(age=age, L0=L0, a=a, b=b, c=c, beta=beta, r=r, am=am, age0=age0)
  La = .QBD_mature(age=age, L0=L0, a=a, b=b, c=c, beta=beta, r=r, am=am, age0=age0)
  L  = Lj + (age>am)*(La - Lj)
  return(L)
}

.QBD_pop = function(age, L0, a, b, c, beta, r, a50, a75, age0=0) {
  p = 0.999 # harcoded
  dx = 0.01 # harcoded
  A = a50*((1-p)/p)^(log(a50/a75)/log(3))
  am = seq(0, A, dx)
  dmat = .dmaturity(am, a50, a75)*dx/p

  L = lapply(am, FUN=.QBD, age=age, L0=L0, a=a, b=b, c=c, beta=beta, r=r, age0=age0)
  L = dmat*do.call(rbind, L)
  L = colSums(L)
  return(L)
}


.Linf_QBD = function(a, b, c, beta) {
  Linf = (c/(a^(1-beta))/r)^(1/(b*(1-beta)))
  return(Linf)
}

# fitting model from simulated data

.c_guess = function(Linf, k, t0, a, b, beta, am, age0) {
  age = seq(age0, am, length.out=10)[-1]
  L0 = .VB(age0, Linf=Linf, k=k, t0=t0)
  L  = .VB(age, Linf=Linf, k=k, t0=t0)
  age = age - age0
  c_guess = ((L^(b*(1-beta)) - L0^(b*(1-beta)))*(a^(1-beta)))/((1-beta)*age)
  return(c_guess)
}

.r_guess = function(a, b, c, beta, Linf) {
  r = c/(a^(1-beta))/(Linf^(b*(1-beta)))
  return(r)
}

lnorm2 = function(obs, sim, tiny=1e-2, ...) {
  if(all(!is.finite(sim))) return(Inf)
  obs = log(obs + tiny)
  sim = log(sim + tiny)
  nlogLike = sum((obs-sim)^2, na.rm=TRUE)
  return(nlogLike)
}

.QBD_par = function(Linf, k, t0, a, b, beta, longevity, a50, a75, age0) {

  age = seq(age0, longevity, by=1/12)

  obs = .VB(x=age, Linf=Linf, k=k, t0=t0)

  cg = .c_guess(Linf=Linf, k=k, t0=t0, a=a, b=b, beta=beta, am=a50, age0=age0)
  rg = .r_guess(a=a, b=b, c=cg, beta=beta, Linf=Linf)
  c_ini = median(cg, na.rm=TRUE)
  r_ini = median(rg, na.rm=TRUE)

  par0 = c(c_ini, r_ini)

  # create function to minimize
  .obj = function(par) {
    sim = .QBD_pop(age=age, L0=L0, a=a, b=b, c=par[1], beta=beta, r=par[2],
                   a50=a50, a75=a75, age0=age0)
    ll = lnorm2(obs=obs, sim=sim)
    return(ll)
  }

  lower = floor(0.8*c(min(cg), min(rg)))
  upper = ceiling(1.2*c(max(cg), max(rg)))

  opt = optim(par=par0, fn=.obj, method="L-BFGS-B", lower=lower, upper=upper)

  sim = .QBD_pop(age=age, L0=L0, a=a, b=b, c=opt$par[1], beta=beta, r=opt$par[2],
                 a50=a50, a75=a75, age0=age0)

  # add everything
  out = c(opt, VB=list(obs), QBD=list(sim))
  return(out)
}


# Reproduction ------------------------------------------------------------

.maturity = function(age, a50, a75) {
  den = 1 + (age/a50)^(log(3)/log(a50/a75))
  out = 1/den
  out[is.nan(out)] = 0
  return(out)
}

.dmaturity = function(age, a50, a75) {
  s = log(3)/log(a50/a75)
  den = 1 + (age/a50)^(s)
  out = -den^(-2)*(age/a50)^(s-1)*s/a50
  out[is.nan(out)] = 0
  return(out)
}





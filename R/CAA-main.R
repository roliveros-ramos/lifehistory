
.simple_caa = function(M, K, alpha.tilde, Y, w, Bini, a_mat, freq=1, ...) {

  A = length(w) # weight control age classes

  maturity = 0 + (seq_len(A) >= a_mat*freq)

  forcing = list(catch=Y, maturity=maturity, selectivity=maturity, weight=w)

  par = list(alpha = .alpha(alpha=alpha.tilde/freq, mortality=M/freq, maturity=maturity, weight=w),
             beta = .beta(alpha=alpha.tilde/freq, K=K, mortality=M/freq, maturity=maturity, weight=w),
             mortality=M/freq,
             Bini=Bini)

  sim = .TVTB(par, forcing)

  return(sim)
}


.TVTB = function(par, forcing) {
  # each par is a matrix T x par_dim
  # this is the basic simulator, outside the parameters are created.
  # forcing
  catch = forcing$catch # tonnes
  ndt = length(catch)
  weight = forcing$weight # grams
  maturity = forcing$maturity
  selectivity = forcing$selectivity
  age_classes = ifelse(is.matrix(weight), ncol(weight), length(weight))
  mat_classes = ifelse(is.matrix(maturity), ncol(maturity), length(maturity))
  sel_classes = ifelse(is.matrix(selectivity), ncol(selectivity), length(selectivity))
  if(age_classes!=mat_classes) stop("Weight and maturity inputs are not consistent.")
  if(age_classes!=sel_classes) stop("Weight and selectivity inputs are not consistent.")
  if(!is.matrix(weight))
    weight = matrix(weight, nrow=ndt, ncol=age_classes, byrow=TRUE)
  if(!is.matrix(maturity))
    maturity = matrix(maturity, nrow=ndt, ncol=age_classes, byrow=TRUE)
  if(!is.matrix(selectivity))
    selectivity = matrix(selectivity, nrow=ndt, ncol=age_classes, byrow=TRUE)
  # parameters
  alpha = par$alpha
  if(length(alpha)==1) alpha = rep(alpha, ndt)
  beta  = par$beta
  if(length(beta)==1) beta = rep(beta, ndt)
  mortality = par$mortality
  if(length(mortality)==1) mortality = rep(mortality, age_classes - 1)
  mort_classes = ifelse(is.matrix(mortality), ncol(mortality), length(mortality)) + 1
  if(age_classes!=mort_classes) stop("Mortality parameters are not consistent with age classes.")
  if(!is.matrix(mortality)) mortality = matrix(mortality, nrow=ndt, ncol=age_classes-1, byrow=TRUE)
  Bini  = par$Bini
  if(length(Bini)==1) {
    Bcomp = RS(mortality[1, ])*weight[1, ]
    Bini = Bini * Bcomp/sum(Bcomp)
  }
  # dimensional checking
  if(length(alpha)!=ndt) stop(sprintf("'alpha' vector must be of length %d.", ndt))
  if(length(beta)!=ndt) stop(sprintf("'beta' vector must be of length %d.", ndt))
  if(length(Bini)!=age_classes)
    stop(sprintf("Initial biomass vector must be of length %d.", age_classes))
  if(nrow(weight)!=ndt) stop(sprintf("'weight' matrix must have %d rows.", ndt))
  if(nrow(maturity)!=ndt) stop(sprintf("'maturity' matrix must have %d rows.", ndt))
  if(nrow(selectivity)!=ndt) stop(sprintf("'selectivity' matrix must have %d rows.", ndt))
  if(nrow(mortality)!=ndt) stop(sprintf("'mortality' parameters are not provided for %d time steps.", ndt))
  # initialization
  B = matrix(NA, nrow=ndt, ncol=age_classes)
  B[1, ] = Bini
  C = matrix(NA, nrow=ndt, ncol=age_classes)
  SSB = rep(NA, ndt)
  # dynamics
  for(t in seq_len(ndt)) {
    SSB[t] = sum(maturity[t, ]*B[t, ])
    b0 = weight[t, 1]*recruitment(ssb=SSB[t], alpha=alpha[t], beta=beta[t], model="Ricker")
    tmp0 =c(b0, exp(diff(log(weight[t, ])))*exp(-mortality[t, ])*B[t, -age_classes])
    cprop = tmp0*selectivity[t, ]/sum(tmp0*selectivity[t, ])
    tmp1 = pmax(tmp0 - catch[t]*cprop, 0.01*tmp0)
    yield = tmp0 - tmp1
    C[t, ] = yield
    if(t!=ndt) B[t+1, ] = tmp1
  }
  # outputs
  output = list(biomass=rowSums(B), yield=rowSums(C), B=B, C=C, SSB=SSB)
  return(output)
}


# Population Dynamics -----------------------------------------------------

RS = function(mortality) {
  # Recruit survival, mortality is a vector for each "age" group, given by freq.
  # mortality has A-1 age classes as plus groups are not supported. A zero is
  # added for trivial 'recruits' survival.
  survival = exp(-c(0, cumsum(mortality)))
  return(survival)
}

LTS = function(mortality, maturity) {
  # LTS: Life time spawners
  # Adult survival: maturity is knife-edge, so RS normalized by first positive.
  survival = RS(mortality) * maturity
  survival = survival / max(survival)
  return(sum(survival))
}

SPR0 = function(mortality, maturity, weight) {
  # SPR0 = recruit survival * weight * maturity
  # weight is given, survival is calculated, maturity is calculated (knife-edge)
  spr = sum(RS(mortality) * maturity * weight)
  return(spr)
}

.alpha =  function(alpha, mortality, maturity, weight) {
  out = alpha*LTS(mortality, maturity)/SPR0(mortality, maturity, weight)
  return(out)
}

alpha = function(alpha, mortality, maturity, weight) {
  out = rep(NA, nrow(mortality))
  for(t in seq_along(out)) {
    out[t] = .alpha(alpha=alpha, mortality=mortality[t, ],
                    maturity=maturity, weight=weight[t, ])
  }
  return(out)
}

.beta = function(alpha, K, mortality, maturity, weight) {
  K = 1e6*K # tonnes -> grammes
  spr0 = SPR0(mortality, maturity, weight)
  out = log(alpha * LTS(mortality, maturity)) * (sum(weight*RS(mortality))/K) / spr0 # 1/gramme
  out = 1e6*out # 1/gramme -> 1/tonnes
  return(out)
}

beta = function(alpha, K, mortality, maturity, weight) {
  out = rep(NA, nrow(mortality))
  for(t in seq_along(out))
    out[t] = .beta(alpha=alpha, K=K[t], mortality=mortality[t, ],
                   maturity=maturity, weight=weight[t, ])
  return(out)
}


ricker = function(ssb, alpha, beta, eta=0) {
  R = alpha*ssb*exp(-beta*ssb)*exp(eta)
  return(R)
}

BH = function(ssb, alpha, beta, eta=0) {
  R = alpha*ssb/(1 + beta*ssb)*exp(eta)
  return(R)
}

recruitment = function(ssb, alpha, beta, eta=0, model="Ricker") {
  # for now, only Ricker recruitment can be used
  R = ricker(ssb=ssb, alpha=alpha, beta=beta, eta=eta)
  return(R)
}

# Regime shifts -----------------------------------------------------------

regime_matrix = function(D, d, T, freq) {

  .regime = function(D, d, T, freq) {
    Dt = D + 0.5*d
    s  = 6/d
    t = head(seq(from=0, to=T, by=1/freq), -1)
    out = 1/(1 + exp(-s*(t - Dt)))
    return(out)
  }

  end = cumsum(D) + head(c(0, d), -1) # actual ends of regimes
  out = cbind(mapply(.regime, D=rev(end), d=rev(d), T=T, freq=freq), 1)
  out = cbind(out[, 1], t(diff(t(out))))
  out = out[, rev(seq_len(ncol(out)))]
  return(out)
}



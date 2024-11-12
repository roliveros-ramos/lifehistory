
a50 = 1
a75 = 1.25

a  = 0.006
b  = 3
beta=0.75

k = 1.2
Linf = 20
t0 = 0
L0 = 6
longevity = 3

age0 = .VB_inv(L0, Linf, k, t0)

opt = .QBD_par(Linf, k, t0, a, b, beta, longevity, a50, a75, age0)

betas = seq(from=0.3, to=0.99, by=0.005)
prof = vector("list", length(betas))
pb = txtProgressBar(style=3, max = length(betas))
for(ibeta in seq_along(betas)) {
  setTxtProgressBar(pb, value=ibeta)
  prof[[ibeta]] = .QBD_par(Linf, k, t0, a, b, betas[ibeta], longevity, a50, a75, age0)
}

c = sapply(prof, FUN=function(x) x$par[1])
r = sapply(prof, FUN=function(x) x$par[2])
ll = sapply(prof, FUN="[[", i="value")

# prof = do.call(rbind, prof)
# c = prof[, 1] #sapply(prof, FUN=function(x) x$par[1])
# r = prof[, 2] #sapply(prof, FUN=function(x) x$par[2])
# ll = prof[, 3] #sapply(prof, FUN="[[", i="value")

sim = sapply(prof, FUN="[[", i="QBD")
obs = prof[[1]]$VB
age = prof[[1]]$age

par(mfrow=c(1, 3), oma=c(1,1,1,1), mar=c(4,4,2,1))
plot(betas, ll, main="log-likelihood profile", type="l", xlab="beta")

matplot(age, sim, type="l", col="grey")
points(age, obs, pch=19)
lines(age, sim[, which.min(ll)], col="blue", lwd=2)
lines(age, sim[, which.min((betas-0.75)^2)], col="red", lwd=2)
lines(age, sim[, which.min((betas-0.66)^2)], col="green", lwd=2)
optlab = sprintf("beta=%0.2f (opt)", betas[which.min(ll)])
mtext(c(optlab, "beta=0.75", "beta=0.66"), side=1,
      col=c("blue", "red", "green"), line=-4:-2, adj=0.95)

plot(betas, c, main="c", type="l", xlab="beta", ylim=c(0, 1.1*max(c)))
lines(betas, r, col="blue")

dsim = apply(sim, 2, diff)

# TOMORROW:
# estimar para todas las cohortes c y r, un solo k.
# resultado: (c,r) por cohorte, un solo k. Optimización única?

# estimation with anchovy 1992-2008.

# anchovy: use c values to estimate temperature dependent parameters.
# sardine, jurel, caballa.

x = out50$bycohort
y = .QBD_fit(x, longevity, a, b, a50, a75, age0)

outx = QBD(out25c, a, b, a50, a75, age0)


library(lifehistory)
library(mvtnorm)

p1 = get_prior("Engraulidae", par=c("Linf", "k"))
p2 = get_prior("Engraulis ringens", par=c("Linf", "k"))
p3 = get_prior("Clupeiformes")
p4 = get_prior("Mustelus mustelus")

x = anchr
MA = 5

ng = .start_grid(x=anchr, MA=5)



rank = "species"
p = p1[[rank]]


par(mfrow=c(3,1), mar=c(3,3,1,1), oma=c(1,1,1,1))

w = c(default=1, class=0, order=0, family=0, genus=0, species=0)
rank = "default"
par$prob0 = dmvnorm(x=par[, c("Linf", "k")], mean=p$species$mean, sigma=p$species$covariance)
par$prob1 = dmvnorm(x=par[, c("Linf", "k")], mean=p$genus$mean, sigma=p$genus$covariance)
par$prob2 = dmvnorm(x=par[, c("Linf", "k")], mean=p$family$mean, sigma=p$family$covariance)
par$prob3 = dmvnorm(x=par[, c("Linf", "k")], mean=p$order$mean, sigma=p$order$covariance)
par$prob4 = dmvnorm(x=par[, c("Linf", "k")], mean=p$class$mean, sigma=p$class$covariance)
par$prob5 = dmvnorm(x=par[, c("Linf", "k")], mean=p$default$mean, sigma=p$default$covariance)

par$prob = with(par, prob0^w[6]*prob1^w[5]*prob2^w[4]*prob3^w[3]*prob4^w[2]*prob5^w[1])

# par$prob = dmvnorm(x=par[, c("Linf", "k")], mean=p[[rank]]$mean,
#                    sigma=p[[rank]]$covariance)


image.plot(exp(Linf), exp(k), prob, las=1)
plot(exp(Linf), pLinf, type="l")
plot(exp(k), pk, type="l")

# function to produce density for bivariate normal from mean and covariance
# function to get prior (uniform) for Linf
Cov = p1$genus$covariance
prob = 0.99

Eigen = eigen(Cov)
Major_radius = sqrt(qchisq(prob, df = 2) * Eigen$values[1])
Minor_radius = sqrt(qchisq(prob, df = 2) * Eigen$values[2])
Angle = atan(Eigen$vector[2, 1]/Eigen$vector[1, 1])/(2*pi)*360

shape::plotellipse
shape::getellipse
shape::rotatexy

shape::plotellipse( rx=Major_radius, ry=Minor_radius, mid=Mean, angle=Angle, lcol=lcol, lty=lty, ...)


p = p2


p = get_prior("Engraulis ringens", par=c("Linf", "k"))



# rank = "species"
# prior = default_prior
# prior = p[[rank]]
# par$prob = dmvnorm(x=log(par[, c("Linf", "k")]),
#                    mean=prior$mean,
#                    sigma=prior$covariance)

# par$obj = par$Rn*par$prob
# prob = matrix(par$prob, nrow=length(Linf), ncol=length(k))
# prob = prob/sum(prob)
# Rn = matrix(par$Rn, nrow=length(Linf), ncol=length(k))
# Rn = Rn/sum(Rn)

# obj = Rn*prob^0.05
# obj = obj/sum(obj)



plot(Linf, pLinf)
plot(k, pk)


par(mfrow=c(3,1))
zlim = range(Rn, prob, obj)*c(0,1.1)
image.plot(Linf, k, Rn, zlim=zlim)
image.plot(Linf, k, prob, zlim=zlim)
image.plot(Linf, k, obj, zlim=zlim)

par[c(which.max(as.numeric(Rn)),
      which.max(as.numeric(prob)),
      which.max(as.numeric(obj))), ]

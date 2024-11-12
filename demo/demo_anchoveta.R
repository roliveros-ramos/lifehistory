library(lubridate)

MA = 3
lifespan = 3

dat = read.csv("data/anchovy.csv", check.names = FALSE)
sel = read.csv("data/anchoveta-empirical_selectivity.csv",
               row.names = 1, check.names = FALSE)[1:31, ]

anch = lf_create(x=t(as.matrix(dat[, -1])), date=dmy(dat$date) + days(14),
                 marks=as.numeric(names(dat)[-1]))

# correcting by selectivity
anchr = list()
anchr[[1]] = window(anch, start=c(1960,1), end=c(1972, 1))/sel[, 2]
anchr[[2]] = window(anch, start=c(1972,1), end=c(1985, 1))/sel[, 3]
anchr[[3]] = window(anch, start=c(1985,1), end=c(1991, 1))/sel[, 4]
anchr[[4]] = window(anch, start=c(1991,1), end=c(1999, 1))/sel[, 5]
anchr[[5]] = window(anch, start=c(1999,1), end=c(2013, 1))/sel[, 6]
anchr[[6]] = window(anch, start=c(2013,1), end=c(2021, 1))/sel[, 7]

anchr = do.call(cbind, anchr)

anch60 = window(anchr, start=c(1960,1), end=c(1972, 1))
anch70 = window(anchr, start=c(1972,1), end=c(1985, 1))
anch80 = window(anchr, start=c(1985,1), end=c(1991, 1))
anch90 = window(anchr, start=c(1991,1), end=c(1999, 1))
anch99 = window(anchr, start=c(1999,1), end=c(2013, 1))
anch13 = window(anchr, start=c(2013,1), end=c(2021, 1))

modes_before = lf_modes(anch)
modes_after  = lf_modes(anchr)

par(mfrow=c(2,1), mar=c(3,3,1,1))
plot(modes_before)
plot(modes_after)

for(n in 1:3) {
  print(n)
  dat = lf_analysis(anchr, MA=2*n+1)
  lmax = tapply(dat$mode[dat$max==1], INDEX = year(dat$date)[dat$max==1], FUN = max)
  print(range(lmax))
  tmax = as.numeric(names(lmax))
}

plot(dat$date, dat$mode, cex=sqrt(dat$peak/max(dat$peak)), pch=19)

xdat = dat[dat$max==1, ]
xdat$year = year(xdat$date)

lmax = aggregate(xdat$mode, by=list(year=xdat$year), FUN=max)
names(lmax)[2] = "lmax"

plot(xdat$date, xdat$mode, cex=sqrt(xdat$peak/max(xdat$peak)),
     las=1, ylim=c(0,20), pch=19)

table(dat$date)

hist(dat$mode[dat$max==1], breaks=30)

x = out30e

age = seq(0, to=x$lifespan, length=100)
out = mapply(FUN=.VB, Linf=x$bycohort$Linf, k=x$bycohort$k, t0=0,
       MoreArgs = list(x=age))
main = .VB(x=age, Linf=x$par$Linf, k=x$par$k, t0=0)
matplot(age, out, type="l", lty=1, col="grey", las=1, ylim=c(0,24))
lines(age, main, col="red", lwd=2)



out5 = ELEFAN1(anchr, lifespan=lifespan, control=list(n=10, MA=7, iter=4), birth = c("01-feb", "01-aug"))

out5 = ELEFAN1(anchr, Linf=c(17,23), k=c(0.5, 2), lifespan=lifespan,
                control=list(n=5, MA=7), birth = c("01-feb", "01-aug"))

out10 = ELEFAN1(anchr, Linf=c(17,23), k=c(0.5, 2), lifespan=lifespan,
                control=list(n=10, MA=7), birth = c("01-feb", "01-aug"), method="rsa")

out30b = ELEFAN1(anchr, Linf=c(15,25), k=c(0.3, 3), lifespan=lifespan,
                control=list(n=30, MA=7), birth = c("01-feb", "01-aug"))

out30c = ELEFAN1(anchr, Linf=c(14,30), k=c(0.3, 2.5), lifespan=lifespan,
                 control=list(n=40, MA=7), birth = c("01-feb", "01-aug"))

out30d = ELEFAN1(anchr, Linf=c(10,40), k=c(0, 2), lifespan=lifespan,
                 control=list(n=c(60, 40), MA=7), birth = c("01-feb", "01-aug"))

out30e = ELEFAN1(anchr, Linf=c(10,40), k=c(0, 2), lifespan=lifespan,
                 control=list(n=c(120, 80), MA=7), birth = c("01-feb", "01-aug"))

out200a = ELEFAN1(anchr, Linf=c(15,23), k=c(0, 1.5), lifespan=lifespan,
                 control=list(n=c(160, 300), MA=7), birth = c("01-feb", "01-aug"))

out10a = ELEFAN1(anchr, Linf=c(15,23), k=c(0, 1.5), lifespan=lifespan,
                  control=list(n=c(8, 16), MA=7), birth = c("01-feb", "01-aug"))

out10a = ELEFAN1(anchr, Linf=c(15,23), k=c(0, 1.5), lifespan=lifespan,
                 control=list(n=c(16, 32), MA=7), birth = c("01-feb", "01-aug"))

out10b = ELEFAN1(anch60, Linf=c(15,23), k=c(0, 1.5), lifespan=lifespan,
                 control=list(n=c(16, 32), MA=7), birth = c("01-feb", "01-aug"))

out10c = ELEFAN1(anch90, Linf=c(15,23), k=c(0, 1.5), lifespan=lifespan,
                 control=list(n=c(16, 32), MA=7), birth = c("01-feb", "01-aug"))


saveRDS(out200a, file="out200.rds")

out60 = ELEFAN1(anchr, Linf=c(17,23), k=c(0.5, 2), lifespan=lifespan,
                control=list(n=60, MA=7), birth = c("01-feb", "01-aug"))

a = 0.0067
b = 3
a50 = .VB_inv(12, Linf=20, k=1.2, t0=0)
a75 = a50 + 0.25
age0 = .VB_inv(6, Linf=20, k=1.2, t0=0)
outx = QBD(out60, a, b, a50, a75, age0)

outx = readRDS("anchoveta-QBD.rds")

zlim = c(0, max(out30e$rsa$Rn, out200a$rsa$Rn))
zlim = range(out30e$rsa$Rn, out200a$rsa$Rn)

image.plot(out30e$rsa$Linf, out30e$rsa$k, out30e$rsa$Rn, zlim=zlim)


par(mfrow=c(1,3), oma=c(1,1,1,1))
image.plot(out10a$rsa$Linf, out10a$rsa$k, out10a$rsa$Rn, zlim=zlim)
image.plot(out30e$rsa$Linf, out30e$rsa$k, out30e$rsa$Rn,
           zlim=zlim, xlim=c(15, 23), ylim=c(0,1.5))
image.plot(out200a$rsa$Linf, out200a$rsa$k, out200a$rsa$Rn, zlim=zlim)

zlim = range(out10a$rsa$Rn, out10b$rsa$Rn, out10c$rsa$Rn)
par(mfrow=c(1,3), oma=c(1,1,1,1))
image.plot(out10a$rsa$Linf, out10a$rsa$k, out10a$rsa$Rn, zlim=zlim)
image.plot(out10b$rsa$Linf, out10b$rsa$k, out10b$rsa$Rn, zlim=zlim)
image.plot(out10c$rsa$Linf, out10c$rsa$k, out10c$rsa$Rn, zlim=zlim)




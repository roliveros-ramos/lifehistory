library(lubridate)
library(fields)
library(parallel)
library(doSNOW)


# TO DO:
# 1. optimize ESP
# 2. optim="", use calibrar.
# 3. time-varying VB.


# 1: given a vector of dates, create the expected length for every cohort present.
# 2: create length for every day, then subset
# 3: melt, filter

dat = read.csv("data/anchovy.csv", check.names = FALSE)

anch = lf_create(x=t(as.matrix(dat[, -1])), date=dmy(dat$date),
            marks=as.numeric(names(dat)[-1]))

# start =
# window(anch$date, start=c(1995,1))

# res = ELEFAN1(x=anch, birth=c("1-feb", "1-ago"),
#               Linf=c(17,23), k=c(0.5, 2), lifespan=lifespan,
#               control=list(n=10, ncores=10))


x = synLFQ5$catch

marks = synLFQ5$midLengths
dates = synLFQ5$dates
lifespan = 7
# obs: k=0.5, Linf=80
MA = 11
x = lf_create(x=x, marks=marks, date=dates)

out10 = ELEFAN1(x, Linf=c(70,90), k=c(0.1, 2), lifespan=lifespan,
                control=list(n=10, MA=MA))

out25 = ELEFAN1(x, Linf=c(60, 110), k=c(0, 4), lifespan=lifespan,
                control=list(n=25, MA=MA))

out25l = ELEFAN1(x, Linf=c(60, 110), k=c(0, 4), lifespan=15,
                control=list(n=25, MA=MA))
out25b = ELEFAN1(x, Linf=c(60, 110), k=c(0.1, 1), lifespan=lifespan,
                control=list(n=25, MA=MA))
out25c = ELEFAN1(x, Linf=c(75, 95), k=c(0.33, 0.66), lifespan=lifespan,
                control=list(n=25, MA=MA))

out25n = ELEFAN1(x, Linf=c(60, 110), k=c(0.1, 1), lifespan=lifespan,
                 control=list(n=50, MA=MA, flag.method="noflag"))
out25m = ELEFAN1(x, Linf=c(60, 110), k=c(0.1, 1), lifespan=lifespan,
                 control=list(n=50, MA=MA, flag.method="mean"))
out25f = ELEFAN1(x, Linf=c(60, 110), k=c(0.1, 1), lifespan=lifespan,
                 control=list(n=50, MA=MA, flag.method="first"))
out25x = ELEFAN1(x, Linf=c(60, 110), k=c(0.1, 1), lifespan=lifespan,
                 control=list(n=50, MA=MA, flag.method="max"))

out50 = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 2), lifespan=lifespan,
                control=list(n=50, MA=MA))
out100 = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 2), lifespan=lifespan,
                 control=list(n=100, MA=MA))
out150 = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 2), lifespan=lifespan,
                 control=list(n=150, MA=MA))
out500 = ELEFAN1(x, Linf=c(60, 100), k=c(0.1, 0.9), lifespan=lifespan,
                 control=list(n=500, MA=MA))
saveRDS(out500, file="out500b.rds")

out500a = readRDS("out500.rds")

out = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 4), lifespan=lifespan, control=list(n=10))

out15 = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 3), lifespan=lifespan, control=list(n=15))
out25 = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 3), lifespan=lifespan, control=list(n=25))

out100 = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 4), lifespan=lifespan, n=100)
out0 = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 4), lifespan=lifespan, n=25, control=control)
out1 = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 4), lifespan=lifespan, n=25,
               control=list(before=TRUE, trim=TRUE, n=25))
plot(out)


plot(out25)
plot(out25b, add=TRUE)
plot(out25c, add=TRUE)

plot(out25b)
plot(out25c)


n = nrow(out25c$bycohort)
age = seq(0, lifespan, by=0.1)
out = NULL
for(i in seq_len(n)) {
  iout = .VB(x=age, k=out25c$bycohort$k[i],
             Linf=out25c$bycohort$Linf[i],
             t0 = out25c$bycohort$t0[i])
  out = cbind(out, iout)
}

image(out25b$rsa$Linf, out25b$rsa$k, out25b$rsa$Rn, add=TRUE)

thr = 0.95
cc = rowSums(x)
cc = cumsum(cc)/sum(cc)
l_min = marks[which.min((cc-thr)^2)]
l_max = 1.2*max(marks)
Linf = range(pretty(c(l_min, l_max)))

cl = makeCluster(5)
registerDoParallel(cl)
out10p = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 3), lifespan=lifespan, control=list(n=10))
# out25p = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 3), lifespan=lifespan, control=list(n=25))
# out15p = ELEFAN1(x, Linf=c(70, 90), k=c(0.1, 3), lifespan=lifespan, control=list(n=15))
stopCluster(cl)


obs
obj = function(t0) {
  sim = cbind(Lt, length=.VB(Lt$age, Linf=Linf, k=k, t0=t0))
  return(-ESP(obs, sim, bins, marks))
}

res = optimize(f=obj, interval = c(-lifespan, 0))





out_a = ELEFAN1(anch90, lifespan=lifespan, control=list(n=30, MA=5),
                birth = c("01-feb", "01-aug"))


out_a2 = ELEFAN1(anch90, lifespan=lifespan, control=list(n=30, MA=5),
                 birth = c("01-feb", "01-aug"))


out_b = ELEFAN1(anch90, lifespan=lifespan, control=list(n=30, MA=5),
                birth = c("01-feb", "01-aug"))

out_c = ELEFAN1(anch90, lifespan=lifespan, control=list(n=30, MA=5),
                birth = c("01-feb", "01-aug"))

beepr::beep(4)

zlim = c(0, max(out_a$Rn, out_b$Rn, out_c$Rn))
col = rev(rainbow(64, start = 0/6, end = 4/6))
col = fields::tim.colors(64)
par(mfrow=c(1,3), oma=c(1,1,1,1))
plot(out_a, zlim=zlim, col=col)
plot(out_b, zlim=zlim, col=col)
plot(out_c, zlim=zlim, col=col)

out = list(out_a$rsa, out_b$rsa, out_c$rsa)

xout = list(rsa=.merge_rsa(out))

zlim = range(xout$rsa$Rn, na.rm=TRUE)
plot.elefan(xout, zlim=zlim, what="Rn")



n = 30
x = out_a2$rsa

k = x$k
Linf = x$Linf

pk = apply(x$Rn, 2, mean)
pLinf = apply(x$Rn, 1, mean)



nk = .new_par(pk, k, 100)
nLinf = .new_par(pLinf, Linf, 100)

plot(k, pk, type="h", col="red")
lines(nk$par, nk$p, type="h", col="blue")

plot(Linf, pLinf, type="h", col="red")
lines(nLinf$par, nLinf$p, type="h", col="blue")


par(mfrow=c(3,1), mar=c(3,3,1,1), oma=c(1,1,1,1))
plot(out_a, zlim=range(0, out_a$rsa$Rn, na.rm=TRUE))
plot(k, log(pk), type="l")
plot(Linf, log(pLinf), type="l")

plot(x$k, pk, type="l")
plot(x$Linf, pLinf, type="l")






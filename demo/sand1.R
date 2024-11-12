library(fields)

data("alba")

rcounts  = lf_restructure(x=alba$catch, MA=5, control=list(trim=TRUE))
rx = lfqRestructure(alba, MA=5)

rr = rcounts/rcounts2
rr[is.na(rr)] = 1

dr = rcounts - rcounts2

par(mfrow=c(1,2))
image.plot(rcounts)
image.plot(rcounts2)

yy = rc$rcounts - rcounts

Fi = c(6, 19, 14, 10, 6, 0, 0, 6, 1, 0, 2, 1)

.lf_restructure_BP1986(Fi, MA=5)
.lf_restructure_GP1997(Fi, MA=5)

i = 6
Fi = alba$catch[, i]

cbind(Fi, 
.lf_restructure_BP1986(Fi, MA=5, trim = FALSE, before=FALSE),
.lf_restructure_BP1986(Fi, MA=5, trim = FALSE, before=TRUE),
.lf_restructure_BP1986(Fi, MA=5, trim = TRUE, before=FALSE),
.lf_restructure_BP1986(Fi, MA=5, trim = TRUE, before=TRUE))

cbind(Fi, 
      .lf_restructure_GP1997(Fi, MA=5, trim = FALSE, before=FALSE),
      .lf_restructure_GP1997(Fi, MA=5, trim = FALSE, before=TRUE),
      .lf_restructure_GP1997(Fi, MA=5, trim = TRUE, before=FALSE),
      .lf_restructure_GP1997(Fi, MA=5, trim = TRUE, before=TRUE))

.lf_restructure_BP1986(Fi, MA=5, trim = FALSE, before=TRUE, lowfreq = TRUE)
.lf_restructure_BP1986(Fi, MA=5, trim = TRUE, before=TRUE, lowfreq = TRUE)
.lf_restructure_BP1986(Fi, MA=5, trim = TRUE, before=FALSE, lowfreq = TRUE)

.lf_restructure_GP1997(Fi, MA=5, trim = FALSE, before=FALSE, lowfreq = TRUE)


x = t(t(Fi))
rcounts  = c(lf_restructure(x, MA=5))
c(lf_restructure(x, MA=5, control=list(trim=TRUE)))
rcounts2 = c(lfqRestructure(list(catch=x), MA=5)$rcounts)

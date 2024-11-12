
library(microbenchmark)

microbenchmark(
ELEFAN(
  lfq = alba,  MA = 7,
  Linf_range = seq(7, 20, length.out = 30),
  K_range = exp(seq(log(0.1),log(4), length.out = 30)),
  method = "cross",
  cross.date = alba$dates[3],
  cross.midLength = alba$midLengths[5],
  contour = TRUE, add.values = FALSE,
  hide.progressbar = TRUE
), 
ELEFAN3(
  lfq = alba,  MA = 7,
  Linf_range = seq(7, 20, length.out = 30),
  K_range = exp(seq(log(0.1),log(4), length.out = 30)),
  method = "cross",
  cross.date = alba$dates[3],
  cross.midLength = alba$midLengths[5],
  contour = TRUE, add.values = FALSE,
  hide.progressbar = TRUE
), times=10)

N = 250

a = ELEFAN(
  lfq = alba,  MA = 7,
  Linf_range = seq(7, 20, length.out = N),
  K_range = exp(seq(log(0.1),log(4), length.out = N)),
  method = "cross", flagging.out = TRUE, 
  cross.date = alba$dates[3],
  cross.midLength = alba$midLengths[5],
  contour = TRUE, add.values = FALSE,
  hide.progressbar = TRUE
)

a3 = ELEFAN3(
  lfq = alba,  MA = 7,
  Linf_range = seq(7, 20, length.out = N),
  K_range = exp(seq(log(0.1),log(4), length.out = N)),
  method = "cross", flagging.out = TRUE, 
  cross.date = alba$dates[3],
  cross.midLength = alba$midLengths[5],
  contour = TRUE, add.values = FALSE,
  hide.progressbar = TRUE
)


score = cbind(ELEFAN=as.numeric(a$score_mat), ELEFAN3=as.numeric(a3$score_mat))
score = as.data.frame(score)
score$diff = score$ELEFAN - score$ELEFAN3
grd = expand.grid(dimnames(a$score_mat))
score$k = as.numeric(as.character(grd$Var1))
score$Linf = as.numeric(as.character(grd$Var2))

score[291,]

alba = lfqRestructure(alba, MA=5)

i = 291
par = list(Linf = score$Linf[i], K = score$k[i], t_anchor = 0.9, C = 0, ts = 0)
x0 = lfqFitCurves(alba, par = par)
x3 = lfqFitCurves2(alba, par = par)
xb = lfqFitCurvesb(alba, par = par)
print(c(x0$ESP, x3$ESP, xb$ESP))


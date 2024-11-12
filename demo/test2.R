library(TropFishR)
library(lubridate)

dat = read.csv("data/anchovy.csv", check.names = FALSE)
ind = rowSums(dat[, -1]) != 0

anchovy = list(sample.no  = seq_len(ncol(dat)-1),
               midLengths = as.numeric(colnames(dat)[-1]),
               dates      = dmy(dat$date)[ind],
               catch      = t(as.matrix(dat[ind, -1])))
class(anchovy) = "lfq"


binsize = unique(diff(anchovy$midLengths))
MA = round(0.23*max(anchovy$midLengths)^0.6/binsize, 0)

kali:::DateStamp()
anchovyt = ELEFAN(lfq = anchovy,  MA = MA,
                  Linf_fix = 20,
  # Linf_range = seq(15, 25, length.out = 10),
  K_range    = exp(seq(log(0.1), log(3), length.out = 10)),
  method     = "optimise",
  contour    = TRUE, add.values = FALSE,
  hide.progressbar = FALSE)
kali:::DateStamp()
plot(anchovyt)

kali:::DateStamp()
anchovyt = ELEFAN(lfq = anchovy,  MA = MA,
                  Linf_range = seq(17, 22, length.out = 100),
                  K_range    = exp(seq(log(0.1), log(3), length.out = 100)),
                  method     = "optimise",
                  contour    = TRUE, add.values = FALSE,
                  hide.progressbar = FALSE)
kali:::DateStamp()
plot(anchovyt)





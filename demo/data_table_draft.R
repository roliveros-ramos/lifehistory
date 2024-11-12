obs = as.data.table(obs)

# brute-force estimation squeleton

kali:::DateStamp()
if(length(Linf)==1) n[1] = 1
Linfs = seq(from=Linf[1], to=Linf[2], length.out=n[1])
ks    = seq(from=k[1], to=k[2], length.out=n[2])
t0s   = seq(from=-1.5, to=0, length.out=150)
par = expand.grid(Linf=Linfs, k=ks, t0=t0s)
sim = cbind(par,
            date=rep(Lt$date, nrow(par)),
            cohort=rep(Lt$cohort, nrow(par)),
            age=rep(Lt$age, nrow(par)))
sim = as.data.table(sim)
sim[, length := .VB(age, Linf, k, t0)]
sim[, length := marks[cut(length, breaks=bins, labels=FALSE)]]
#
sim = obs[sim, on=c("date", "length")]
sim = sim[!is.na(length)]
# setkey(sim, ID)
neg = sim[cohort_peaks==0, .(neg=sum(Fs, na.rm=TRUE)), by=.(Linf, k, t0)]
pos = sim[cohort_peaks!=0, .SD[which.min(cohort), .(Fs)],
          by=.(Linf, k, t0, cohort_peaks), .SDcols=c("cohort", "Fs")]
pos = pos[, .(pos=sum(Fs)), by=.(Linf, k, t0)]

par = as.data.table(par)

esp = neg[pos[par, on=c("Linf", "k", "t0")], , on=c("Linf", "k", "t0")]
esp[is.na(pos), pos := 0]
esp[is.na(neg), neg := 0]
esp = esp[, .(Linf, k, t0, esp=neg+pos)]

rsa = esp[, .SD[which.max(esp), .(esp)], by=.(Linf, k)]

mm = matrix(rsa$esp, nrow=length(Linfs), ncol=length(ks), byrow = TRUE)
mm = 0.1*10^(mm/peaks$ASP)

xout = list(rsa=list(Rn=mm, Linf=Linfs, k=ks))
class(xout) = "elefan"
kali:::DateStamp()
plot(xout)

image.plot(Linfs, ks, mm)

?data.table

esp = ESP2(obs, sim, bins, marks)



# test --------------------------------------------------------------------


if(length(Linf)==1) n[1] = 1
Linfs = seq(from=Linf[1], to=Linf[2], length.out=n[1])
ks    = seq(from=k[1], to=k[2], length.out=n[2])
t0s   = seq(from=-1, to=0, length.out=100)
par = expand.grid(Linf=Linfs, k=ks, t0=t0s)
sim = data.table(ID=rep(as.integer(seq_len(nrow(par))), each=nrow(Lt)),
                 date=rep(Lt$date, nrow(par)),
                 cohort=rep(Lt$cohort, nrow(par)),
                 age=rep(Lt$age, nrow(par)))
sim[, length := .VB(age, par$Linf, par$k, par$t0)]
sim[, length := marks[cut(length, breaks=bins, labels=FALSE)]]
#
sim = obs[sim, on=c("date", "length")]
sim = sim[!is.na(length)]
# setkey(sim, ID)
neg = sim[cohort_peaks==0, .(neg=sum(Fs, na.rm=TRUE)), by=ID]
kali:::DateStamp()
pos = sim[cohort_peaks!=0, .SD[which.min(cohort), .(Fs)],
          by=.(ID, cohort_peaks), .SDcols=c("cohort", "Fs")]
kali:::DateStamp()
pos = pos[, .(pos=sum(Fs, na.rm=TRUE)), by=ID]

esp = pos[neg, on="ID"]
esp[is.na(pos), pos := 0]
setkey(esp, "ID")


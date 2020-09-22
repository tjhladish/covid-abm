d = read.csv("plot_log.csv")
d$date = as.Date(d$date)
#ed = read.csv("rcasedeath-florida.csv")
#ed = read.csv("rcasedeath-escambia.csv")
ed = read.csv("rcasedeath-dade.csv")
ed$Date = as.Date(ed$Date)

ticks <- seq(as.Date('2020-03-01'), d$date[length(d$date)]+5, by = "months") # so much stupid

png('simvis.png', height=1200, width=2400, res=200)
par(mfrow=c(4,2), mar=c(1,4.1,0.5,1), oma=c(2,0,1,0))

plot(d$date, d$sd, col='darkorange3', ylim=c(0,1), type='l', xlab='', ylab='Interventions', xaxt='n')
abline(v=d$date[d$closed==1], col='#f0e68c55', lwd=5)
axis(1, at=ticks, labels=F)

plot(d$date, d$cinf, col='darkslateblue', type='l', xlab='', ylab='Cumulative infections', xaxt='n')
axis(1, at=ticks, labels=F)

plot(d$date, d$Rt, col='red3', type='l', xlab='', ylab=expression(R[t]), xaxt='n')
abline(h=1.0, lty=3)
axis(1, at=ticks, labels=F)

d$crcase   = cumsum(d$rcase)
d$crdeath  = cumsum(d$rdeath)
ed$rcase = ed$rcase
ed$rdeath = ed$rdeath
ed$crcase  = cumsum(ed$rcase)
ed$crdeath = cumsum(ed$rdeath)

plot(d$date, d$crcase/max(d$crcase), col='royalblue3', type='l', xlab='', ylab='Sim CRC & CRD', xaxt='n', ylim=c(0,1))
lines(d$date, d$crdeath/max(d$crdeath), col='green4')

ymax = max(d$rcase, ed$rcase)
plot(d$date, d$rcase, col='royalblue3', type='l', xlab='', ylab='Reported cases', xaxt='n', ylim=c(0,ymax))
lines(ed$Date, ed$rcase)
axis(1, at=ticks, labels=F)

ymax = max(d$crcase, ed$crcase)
plot(d$date, d$crcase, col='royalblue3', type='l', xlab='', ylab='Reported cases', xaxt='n', ylim=c(0,ymax))
lines(ed$Date, ed$crcase)
axis(1, at=ticks, labels=F)

ymax = max(d$rdeath, ed$rdeath)
plot(d$date, d$rdeath, col='green4', type='l', xlab='', ylab='Reported deaths', ylim=c(0,ymax))
lines(ed$Date, ed$rdeath)

ymax = max(d$crdeath, ed$crdeath)
plot(d$date, d$crdeath, col='green4', type='l', xlab='', ylab='Reported deaths', ylim=c(0,ymax))
lines(ed$Date, ed$crdeath)

#plot(d$date, d$rcase, col='royalblue3', type='l', xlab='')
#lines(ed$Date), ed$rcase)
#
#plot(d$date, d$rdeath, col='green4', type='l', xlab='')
#lines(ed$Date), ed$rdeath)
dev.off()

d = read.csv("plot_log.csv")
d$date = as.Date(d$date)
ed = read.csv("tschart-escambia.csv")
ed$ChartDate = as.Date(ed$ChartDate)
#ed = read.csv("ts-escambia.csv"); ed$ChartDate = ed$EventDate 

ticks <- seq(as.Date('2020-03-01'), d$date[length(d$date)]+5, by = "months") # so much stupid

png('simvis.png', height=1200, width=1800, res=200)
par(mfrow=c(5,1), mar=c(1,4.1,0.5,1), oma=c(2,0,1,0))

plot(d$date, d$sd, col='darkorange3', ylim=c(0,1), type='l', xlab='', ylab='Interventions', xaxt='n')
abline(v=d$date[d$closed==1], col='#f0e68c55', lwd=5)
axis(1, at=ticks, labels=F)

plot(d$date, d$Rt, col='red3', type='l', xlab='', ylab=expression(R[t]), xaxt='n')
abline(h=1.0, lty=3)
axis(1, at=ticks, labels=F)

plot(d$date, d$cinf, col='darkslateblue', type='l', xlab='', ylab='Cumulative infections', xaxt='n')
axis(1, at=ticks, labels=F)

plot(d$date, cumsum(d$rcase), col='royalblue3', type='l', xlab='', ylab='Reported cases', xaxt='n')
lines(ed$ChartDate, cumsum(ed$rcase))
axis(1, at=ticks, labels=F)

plot(d$date, cumsum(d$rdeath), col='green4', type='l', xlab='', ylab='Reported deaths')
lines(ed$ChartDate, cumsum(ed$rdeath))

#plot(d$date, d$rcase, col='royalblue3', type='l', xlab='')
#lines(ed$ChartDate), ed$rcase)
#
#plot(d$date, d$rdeath, col='green4', type='l', xlab='')
#lines(ed$ChartDate), ed$rdeath)
dev.off()

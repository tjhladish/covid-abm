calc_centered_avg = function(vals, halfwindow) {
    n       = length(vals)
    indices = 1:n
    ends    = pmin(indices + halfwindow + 1, n)
    starts  = ifelse(halfwindow >= indices, 0, indices - halfwindow)
    widths  = ends - starts
    ma      = mapply(function(s,e) sum(vals[s:e]), s = starts, e = ends) / widths
    return(ma)
}

d = read.csv("plot_log.csv")
d$date = as.Date(d$date)
#escambia_fraction    = 0.0153 # fraction of FL pop that lives in Escambia
pop_florida  = 20609673
pop_escambia = 312212
pop_dade     = 2794464

ed = read.csv("rcasedeath-florida.csv"); per10k = 1e4/pop_florida
#ed = read.csv("rcasedeath-escambia.csv"); per10k = 1e4/pop_escambia
#ed = read.csv("rcasedeath-dade.csv"); per10k = 1e4/pop_dade
ed$Date = as.Date(ed$Date)

ticks <- seq(as.Date('2020-03-01'), d$date[length(d$date)]+5, by = "months") # so much stupid


#        {"name" : "mar2apr30_rcases",  "short_name" : "rcases1",    "num_type" : "FLOAT", "value" :  16.7},
#        {"name" : "jun1oct5_rcases",   "short_name" : "rcases2",    "num_type" : "FLOAT", "value" :  322},
#        {"name" : "mar2oct5_rcases",   "short_name" : "rcases_all", "num_type" : "FLOAT", "value" :  349},
#        {"name" : "jun15jul05_slope",  "short_name" : "slope1",     "num_type" : "FLOAT", "value" :  0.160},
#        {"name" : "jul25aug14_slope",  "short_name" : "slope2",     "num_type" : "FLOAT", "value" :  -0.106},
#        {"name" : "mar2may31_rdeaths", "short_name" : "rdeaths1",   "num_type" : "FLOAT", "value" :  1.23},
#        {"name" : "jun1oct5_rdeaths",  "short_name" : "rdeath2",    "num_type" : "FLOAT", "value" :  6.02},
#        {"name" : "mar15sep5_deaths",  "short_name" : "deaths",     "num_type" : "FLOAT", "value" :  9.75},
#        {"name" : "mar18oct5_hosps",   "short_name" : "hosps",      "num_type" : "FLOAT", "value" :  21.7},


png('simvis.png', height=1200, width=2400, res=240)
par(mfrow=c(4,2), mar=c(1,4.1,0.5,1), oma=c(2,0,1,0))

# Transmission reductions
plot(d$date, d$sd, col='darkorange3', ylim=c(0,1), type='l', xlab='', ylab='Interventions', xaxt='n')
abline(v=d$date[d$closed==1], col='#f0e68c55', lwd=5)
axis(1, at=ticks, labels=F)

# Cumulative infections
plot(d$date, d$cinf, col='darkslateblue', type='l', xlab='', ylab='Cumulative infections', xaxt='n')
axis(1, at=ticks, labels=F)
abline(h=d$cinf[d$date == as.Date('2020-06-01')], lty=2)
abline(v=as.Date('2020-06-01'), lty=2)

# Rt
plot(d$date, d$Rt, col='red3', type='l', xlab='', ylab=expression(R[t]), xaxt='n')
points(as.Date(c('2020-03-27', '2020-05-07', '2020-06-28')), c(0.98, 1.0, 1.0))
segments(as.Date(c('2020-03-27', '2020-05-07')), c(0.78, 1.37), as.Date(c('2020-05-07', '2020-06-28')))
abline(h=1.0, lty=3)
axis(1, at=ticks, labels=F)

d$crcase   = cumsum(d$rcase)
d$crdeath  = cumsum(d$rdeath)
#death_underreporting = 20.1/11.8 # excess death / recognized COVID deaths
ed$rcase = ed$rcase*per10k
ed$rdeath = ed$rdeath*per10k #*death_underreporting
ed$crcase  = cumsum(ed$rcase)
ed$crdeath = cumsum(ed$rdeath)

# probably get rid of this
plot(d$date, d$crcase/max(d$crcase), col='royalblue3', type='l', xlab='', ylab='Sim CRC & CRD', xaxt='n', ylim=c(0,1))
lines(d$date, d$crdeath/max(d$crdeath), col='green4')

#        {"name" : "jul19_rcases",      "short_name" : "jul19",      "num_type" : "FLOAT", "value" :  5.26},
#        {"name" : "peak_julian_day",   "short_name" : "peak_day",   "num_type" : "INT",   "value" :  201},
# Reported cases
ymax = max(d$rcase, ed$rcase)
plot(d$date, d$rcase, col='royalblue3', type='l', xlab='', ylab='Reported cases', xaxt='n', ylim=c(0,ymax))
lines(ed$Date, ed$rcase)
axis(1, at=ticks, labels=F)

ymax = max(d$crcase, ed$crcase)
plot(d$date, d$crcase, col='royalblue3', type='l', xlab='', ylab='Reported cases', xaxt='n', ylim=c(0,ymax))
lines(ed$Date, ed$crcase)
axis(1, at=ticks, labels=F)

ymax = max(d$rdeath, ed$rdeath)
plot(d$date, d$rdeath, col='green4', type='l', xlab='', ylab='Reported deaths', xaxt='n', ylim=c(0,ymax))
lines(ed$Date, ed$rdeath)
axis(1, at=ticks, labels=format(ticks, "%b"))

ymax = max(d$crdeath, ed$crdeath)
plot(d$date, d$crdeath, col='green4', type='l', xlab='', ylab='Reported deaths', xaxt='n', ylim=c(0,ymax))
lines(ed$Date, ed$crdeath)
axis(1, at=ticks, labels=format(ticks, "%b"))

#plot(d$date, d$rcase, col='royalblue3', type='l', xlab='')
#lines(ed$Date), ed$rcase)
#
#plot(d$date, d$rdeath, col='green4', type='l', xlab='')
#lines(ed$Date), ed$rdeath)
dev.off()

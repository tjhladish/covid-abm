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
ed$date = as.Date(ed$Date)

ticks <- seq(as.Date('2020-03-01'), d$date[length(d$date)]+5, by = "months") # so much stupid

# for shading the background, we need an even number of months; add 1 if necessary
shaded_boundaries = ticks
if (length(shaded_boundaries) %% 2 == 1) { shaded_boundaries = c(shaded_boundaries, seq(shaded_boundaries[length(shaded_boundaries)], length=2, by="1 month")[2]); }

shading = function() {
    shaded_colors = c('#EEEEEE', '#EEFFEE', '#EEEEFF')
    rect(shaded_boundaries[c(T,F)], par("usr")[3], shaded_boundaries[c(F,T)], par("usr")[4], col=shaded_colors, lwd=0)
}

annotate = function(text) {
    legend('topleft', text, bty='n')
    box()
    axis(1, at=ticks, labels=F)
    axis(1, at=ticks+15, tick=F, padj=-1.5, gap.axis=0.25, labels=substring(format(ticks, "%b"), 1, 1))
#axis(1, at=ticks, labels=format(ticks, "%b"))
}

#        {"name" : "mar2apr30_rcases",  "short_name" : "rcases1",    "num_type" : "FLOAT", "value" :  16.7},
#        {"name" : "jun1oct5_rcases",   "short_name" : "rcases2",    "num_type" : "FLOAT", "value" :  322},
#        {"name" : "mar2oct5_rcases",   "short_name" : "rcases_all", "num_type" : "FLOAT", "value" :  349},
#        {"name" : "jun15jul05_slope",  "short_name" : "slope1",     "num_type" : "FLOAT", "value" :  0.160},
#        {"name" : "jul25aug14_slope",  "short_name" : "slope2",     "num_type" : "FLOAT", "value" :  -0.106},
#        {"name" : "mar2may31_rdeaths", "short_name" : "rdeaths1",   "num_type" : "FLOAT", "value" :  1.23},
#        {"name" : "jun1oct5_rdeaths",  "short_name" : "rdeath2",    "num_type" : "FLOAT", "value" :  6.02},
#        {"name" : "mar15sep5_deaths",  "short_name" : "deaths",     "num_type" : "FLOAT", "value" :  9.75},
#        {"name" : "mar18oct5_hosps",   "short_name" : "hosps",      "num_type" : "FLOAT", "value" :  21.7},


png('simvis.png', height=1200, width=2400, res=200)
par(mfrow=c(3,3), mar=c(1,2.1,0.5,0.5), oma=c(2,0.5,1,0.5))

# Social distancing
plot(d$date, d$sd, col='darkorange3', ylim=c(0,1), type='n', xlab='', ylab='', xaxt='n', bty='n')
shading()
abline(v=d$date[d$closed==1], col='#f0e68c55', lwd=5)
lines(d$date, d$sd, col='darkorange3')
annotate('Social distancing')

# Rt
plot(d$date, d$Rt, type='n', xlab='', ylab='', xaxt='n', bty='n')
shading()
lines(d$date, d$Rt, col='red3')
points(as.Date(c('2020-03-27', '2020-05-07', '2020-06-28')), c(0.98, 1.0, 1.0))
segments(as.Date(c('2020-03-27', '2020-05-07')), c(0.78, 1.37), as.Date(c('2020-05-07', '2020-06-28')))
abline(h=1.0, lty=3)
annotate(expression(R[t]))

# Cumulative infections
plot(d$date, d$cinf, type='n', xlab='', ylab='', xaxt='n', bty='n')
shading()
lines(d$date, d$cinf, col='darkslateblue')
annotate('Cumulative infections')

# Seasonality
plot(d$date, d$seasonality, type='n', xlab='', ylab='', xaxt='n', bty='n')
shading()
abline(h=1, lty=3)
lines(d$date, d$seasonality, col='purple')
annotate('Seasonality')

d$crcase   = cumsum(d$rcase)
d$crdeath  = cumsum(d$rdeath)
#death_underreporting = 20.1/11.8 # excess death / recognized COVID deaths
ed$rcase = ed$rcase*per10k
#ed$rhosp = ed$rhosp*per10k
ed$rdeath = ed$rdeath*per10k #*death_underreporting
ed$crcase  = cumsum(ed$rcase)
ed$crdeath = cumsum(ed$rdeath)

#plot(d$date, d$crcase/max(d$crcase), col='royalblue3', type='l', xlab='', ylab='Sim CRC & CRD', xaxt='n', ylim=c(0,1))
#lines(d$date, d$crdeath/max(d$crdeath), col='green4')

# VOC prevalence
plot(d$date, d$vocprev1, type='n', xlab='', ylab='', xaxt='n', bty='n')
shading()
lines(d$date, d$vocprev1, col='royalblue3')
lines(d$date, d$vocprev2, col='turquoise4')
annotate('VOC prevalence')

#        {"name" : "jul19_rcases",      "short_name" : "jul19",      "num_type" : "FLOAT", "value" :  5.26},
#        {"name" : "peak_julian_day",   "short_name" : "peak_day",   "num_type" : "INT",   "value" :  201},

# cumulative reported cases
ymax = max(d$crcase, ed$crcase, na.rm=T)
plot(d$date, d$crcase, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
shading()
lines(ed$date, ed$crcase)
lines(d$date, d$crcase, col='royalblue3')
annotate('Cumulative reported cases')

# Reported cases
ymax = max(d$rcase, ed$rcase, na.rm=T)
plot(d$date, d$rcase, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
shading()
lines(ed$date, ed$rcase)
lines(d$date, d$rcase, col='royalblue3')
annotate('Reported cases')

# reported deaths
ymax = max(d$rdeath, ed$rdeath, na.rm=T)
plot(d$date, d$rdeath, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
shading()
lines(ed$date, ed$rdeath)
lines(d$date, d$rdeath, col='green4')
annotate('Reported deaths')
#axis(1, at=ticks, labels=format(ticks, "%b"))

# cumulative reported deaths
ymax = max(d$crdeath, ed$crdeath, na.rm=T)
plot(d$date, d$crdeath, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
shading()
lines(ed$date, ed$crdeath)
lines(d$date, d$crdeath, col='green4')
annotate('Cumulative reported deaths')

#plot(d$date, d$rcase, col='royalblue3', type='l', xlab='')
#lines(ed$date), ed$rcase)
#
#plot(d$date, d$rdeath, col='green4', type='l', xlab='')
#lines(ed$date), ed$rdeath)
dev.off()

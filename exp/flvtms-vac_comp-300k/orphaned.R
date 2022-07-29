
















ticks <- seq(as.Date('2020-02-01'), d$date[length(d$date)]+5, by = "months") # so much stupid

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

png(paste0('simvis', suffix, '.png'), height=1200, width=2400, res=200)
par(mfrow=c(5,3), mar=c(1,2.1,0.5,0.5), oma=c(2,0.5,1,0.5))

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
ymax = max(d$cinf, d$seroprev, seroprev$upper, 1, na.rm=T)
plot(d$date, d$cinf, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
shading()
abline(h=1, lty=3)
arrows(seroprev$date, seroprev$lower, seroprev$date, seroprev$upper, length = 0.05, code=0, col = rgb(0,0,0,0.25))
points(seroprev$date, seroprev$est, pch = 20, cex = 0.5)
lines(d$date, d$cinf, col='orangered', lty = 3)
lines(d$date, d$seroprev, col='orangered')
annotate('Seroprevalence (solid) and cumulative infections (dashed)')

# Seasonality
plot(d$date, d$seasonality, type='n', xlab='', ylab='', xaxt='n', bty='n')
shading()
abline(h=1, lty=3)
lines(d$date, d$seasonality, col='purple')
annotate('Seasonality')

#plot(d$date, d$crcase/max(d$crcase), col='royalblue3', type='l', xlab='', ylab='Sim CRC & CRD', xaxt='n', ylim=c(0,1))
#lines(d$date, d$crdeath/max(d$crdeath), col='green4')

# VOC prevalence
plot(d$date, d$vocprev1, type='n', xlab='', ylab='', xaxt='n', bty='n')
shading()
lines(d$date, d$vocprev1, col='royalblue3')
lines(d$date, d$vocprev2, col='turquoise4')
lines(d$date, d$vocprev3, col='darkorchid3')
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
abline(h=0)
lines(ed$date, ed$rcase)
lines(d$date, d$rcase, col='royalblue3')
annotate('Reported cases')

# reported deaths
ymax = max(d$rdeath, ed$rdeath, na.rm=T)
plot(d$date, d$rdeath, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
shading()
lines(ed$date, ed$rdeath)
lines(d$date, d$rdeath, col='green4')
annotate('All deaths')
#axis(1, at=ticks, labels=format(ticks, "%b"))

# cumulative reported deaths
ymax = max(d$crdeath, ed$crdeath, na.rm=T)
plot(d$date, d$crdeath, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
shading()
lines(ed$date, ed$crdeath)
lines(d$date, d$crdeath, col='green4')
annotate('Cumulative (excess) deaths')

# infections by vax status
ymax = max(d$vaxInfs + d$unvaxInfs, na.rm=T)
plot(d$date, d$vaxInfs, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
shading()
lines(d$date, d$vaxInfs, col='dodgerblue4')
lines(d$date, d$unvaxInfs + d$vaxInfs, col='coral')
annotate('Infections by vax status (stacked)')

# VES over time
plot(d$date, d$VES, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,1), bty='n')
shading()
lines(d$date, d$VES, col='dodgerblue1')
lines(d$date, d$VESAvg, col='dodgerblue4', lwd = 2)
annotate('Direct VES')

# Vaccine delivery
ymax = max(vax$cov1, d$cov1, 1, na.rm=T)
plot(d$date, d$std_doses, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
shading()
abline(h=1, lty=3)
lines(vax$date, vax$cov1, lty = 3)
lines(vax$date, vax$cov2)
lines(d$date, d$cov1, col='purple', lty=3)
lines(d$date, d$cov2, col='purple')
annotate('Vaccination coverage (dose 1 - dashed; dose 2 - solid)')



# icu inc/prev
# ymax = max(d$icuInc, d$icuPrev, rm.na=T)
# plot(d$date, d$icuPrev, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
# shading()
# lines(d$date, d$icuPrev, col='orangered')
# lines(d$date, d$icuInc, col='orangered', lty=3)
# annotate('ICU prev (solid) and inc (dotted)')

# breakthrough ratio over time
plot(d$date, d$brkthruRatio, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,1), bty='n')
shading()
lines(cdc$date, cdc$brkthruRatio)
lines(d$date, d$brkthruRatio, col='tan')
lines(d$date, d$brkthruRatioAvg, col='tan4', lwd = 2)
annotate('Breakthrough ratio')

#plot(d$date, d$rcase, col='royalblue3', type='l', xlab='')
#lines(ed$date), ed$rcase)
#
#plot(d$date, d$rdeath, col='green4', type='l', xlab='')
#lines(ed$date), ed$rdeath)
invisible(dev.off())

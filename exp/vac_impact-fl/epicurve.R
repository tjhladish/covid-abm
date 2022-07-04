# needs simvis.R to be called first
# Reported cases
ymax = max(d$rcase, ed$rcase, na.rm=T)
png('projected_cases.png', width=1800, height=1000, res=180)
par(mar=c(5.1, 4.1, 0.5, 0.5))
plot(d$date, d$rcase, type='n', xlab='Date (month)', xaxt='n', ylab='Reported cases per 10k', ylim=c(0,ymax), bty='n')
shading()
lines(ed$date, ed$rcase, lwd=2)
lines(d$date, d$rcase, col='red')
annotate('')
legend('topleft', legend=c('Reported cases', 'Modeled cases'), col=c('black', 'red'), lty=1, lwd=3)
dev.off()

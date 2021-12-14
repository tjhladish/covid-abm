# needs simvis.R to be called first

fl_rescale = 2148 # pop of FL/10k

# CURRENTLY ASSUMES SIM DATA STARTS ON WEDNESDAY, FEB 5, 2020
full_weeks = floor((nrow(d) - 2)/7)
remainder  = (nrow(d) - 2) %% 7
week_idx = c(c(0,0), rep(1:full_weeks, each=7), rep(full_weeks+1, remainder)) #
d$week_idx = week_idx
weekly = aggregate(d[,c('rcase', 'rdeath')], by=list(week_idx=d$week_idx), FUN=sum)
ed_tmp = merge(ed, d[,c('date', 'week_idx')], by='date')
ed_weekly = aggregate(ed_tmp[,c('rcase', 'rdeath')], by=list(week_idx=ed_tmp$week_idx), FUN=sum)
if (sum(ed_tmp$week_idx == min(ed_tmp$week_idx)) < 7) { ed_weekly = ed_weekly[-1,] }
if (sum(ed_tmp$week_idx == max(ed_tmp$week_idx)) < 7) { ed_weekly = ed_weekly[1:(nrow(ed_weekly)-1),] }
weekly$rdate = as.Date('2020-02-07') + (7 * 0:(nrow(weekly) - 1))
names(ed_weekly) = c('week_idx', 'ercase', 'erdeath')
weekly = merge(weekly, ed_weekly, by='week_idx', all.x=T)


# Daily reported cases
ymax = max(d$rcase*fl_rescale, ed$rcase*fl_rescale, na.rm=T)
ymax = 32894.47
png('projected_cases_211006.png', width=1800, height=900, res=180)
par(mar=c(3.1, 4.1, 0.5, 0.5))
plot(d$date, d$rcase*fl_rescale, type='n', xlab='', xaxt='n', ylab='Reported cases per day', ylim=c(0,ymax), bty='n')
title(xlab='Date (month)', line=2)
shading()
lines(ed$date, ed$rcase*fl_rescale, lwd=2)
lines(d$date, d$rcase*fl_rescale, col='orange')
annotate('')
legend('topleft', legend=c('Reported cases', 'Modeled reported cases'), col=c('black', 'orange'), lty=1, lwd=3)
dev.off()


# Weekly reported cases
ymax = max(weekly$rcase*fl_rescale, weekly$ercase*fl_rescale, na.rm=T)
png('projected_cases-weekly_211006.png', width=1800, height=900, res=180)
par(mar=c(3.1, 4.1, 0.5, 0.5))
plot(weekly$rdate, weekly$rcase*fl_rescale, type='n', xlab='', xaxt='n', ylab='Reported cases per week', ylim=c(0,ymax), bty='n')
title(xlab='Date (month)', line=2)
shading()
lines(weekly$rdate, weekly$ercase*fl_rescale, lwd=2)
lines(weekly$rdate, weekly$rcase*fl_rescale, col='orange')
annotate('')
legend('topleft', legend=c('Reported cases', 'Modeled reported cases'), col=c('black', 'orange'), lty=1, lwd=3)
dev.off()


# Rt
png('projected_Rt_211006.png', width=1800, height=900, res=180)
par(mar=c(3.1, 4.1, 0.5, 0.5))
plot(d$date, d$Rt, type='n', xlab='', ylab=expression(R[t]), xaxt='n', bty='n')
title(xlab='Date (month)', line=2)
shading()
lines(d$date[1:(length(d$date)-9)], d$Rt[1:(length(d$date)-9)], col='red3')
#points(as.Date(c('2020-03-27', '2020-05-07', '2020-06-28')), c(0.98, 1.0, 1.0))
#segments(as.Date(c('2020-03-27', '2020-05-07')), c(0.78, 1.37), as.Date(c('2020-05-07', '2020-06-28')))
abline(h=1.0, lty=3)
annotate('')
dev.off()


# All infections 
#ymax = max(d$rcase, ed$rcase, na.rm=T)
png('projected_infections_211006.png', width=1800, height=900, res=180)
par(mar=c(3.1, 4.1, 0.5, 0.5))
plot(d$date, d$inf*fl_rescale, type='n', xlab='', xaxt='n', ylab='Modeled infections per day', bty='n')
title(xlab='Date (month)', line=2)
shading()
lines(d$date[-1], d$inf[-1]*fl_rescale, col='green4')
annotate('')
#legend('topleft', legend=c('', 'Modeled cases'), col=c('black', 'red'), lty=1, lwd=3)
dev.off()


# Reported deaths
ymax = max(d$rdeath*fl_rescale, ed$rdeath*fl_rescale, na.rm=T)
png('projected_deaths_211006.png', width=1800, height=900, res=180)
par(mar=c(3.1, 4.1, 0.5, 0.5))
plot(d$date, d$rdeath*fl_rescale, type='n', xlab='', xaxt='n', ylab='Reported deaths per day', ylim=c(0,ymax), bty='n')
title(xlab='Date (month)', line=2)
shading()
lines(ed$date, ed$rdeath*fl_rescale, lwd=2)
lines(d$date, d$rdeath*fl_rescale, col='red')
annotate('')
legend('topleft', legend=c('Reported deaths', 'Modeled reported deaths'), col=c('black', 'red'), lty=1, lwd=3)
dev.off()


# Reported hospitalizations
#ymax = max(d$rhosp, ed$rhosp, na.rm=T)
#png('projected_hospitalizations_211006.png', width=1800, height=900, res=180)
#par(mar=c(3.1, 4.1, 0.5, 0.5))
#plot(d$date, d$rhosp, type='n', xlab='', xaxt='n', ylab='Reported hospitalizations', ylim=c(0,ymax), bty='n')
#title(xlab='Date (month)', line=2)
#shading()
#lines(ed$date, ed$rhosp, lwd=2)
#lines(d$date, d$rhosp, col='red')
#annotate('')
#legend('topleft', legend=c('Reported hospitalizations', 'Modeled hospitalizations'), col=c('black', 'red'), lty=1, lwd=3)
#dev.off()

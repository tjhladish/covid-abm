#!/usr/bin/env Rscript
library(grDevices)

args = commandArgs(trailingOnly=TRUE)

suffix = '-scenario3'
if (length(args)>0) {suffix = args[1]}

calc_centered_avg = function(vals, halfwindow) {
    n       = length(vals)
    indices = 1:n
    ends    = pmin(indices + halfwindow + 1, n)
    starts  = ifelse(halfwindow >= indices, 0, indices - halfwindow)
    widths  = ends - starts
    ma      = mapply(function(s,e) sum(vals[s:e]), s = starts, e = ends) / widths
    return(ma)
}

d = read.csv(paste0("plot_log", suffix, ".csv"))
d$date = as.Date(d$date)

#escambia_fraction    = 0.0153 # fraction of FL pop that lives in Escambia
pop_florida  = 21538187 #20609673
pop_escambia = 312212
pop_dade     = 2794464

ed = read.csv("rcasedeath-florida.csv"); per10k = 1e4/pop_florida
#ed = read.csv("rcasedeath-escambia.csv"); per10k = 1e4/pop_escambia
#ed = read.csv("rcasedeath-dade.csv"); per10k = 1e4/pop_dade
ed$date = as.Date(ed$Date)

cdcFull = read.csv("Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Vaccination_Status.csv", stringsAsFactors=F)
cdc = cdcFull[cdcFull$Age.group == "all_ages_adj" & cdcFull$Vaccine.product == "all_types" & cdcFull$outcome == "case",]
cdc$date = as.Date("2020-12-27") + (cdc$MMWR.week*7)

d$crcase   = cumsum(d$rcase)
d$crdeath  = cumsum(d$rdeath)
is.na(d) = sapply(d, is.infinite)
d$VESAvg = filter(d$VES, rep(1/7, 7), sides = 2)
d$brkthruRatioAvg = filter(d$brkthruRatio, rep(1/7, 7), sides = 2)
#death_underreporting = 20.1/11.8 # excess death / recognized COVID deaths
ed$rcase = ed$rcase*per10k
#ed$rhosp = ed$rhosp*per10k
#ed$rdeath = ed$rdeath*per10k #*death_underreporting
ed$rdeath = ed$death_incd*per10k
ed$crcase  = cumsum(ed$rcase)
ed$crdeath = cumsum(ed$rdeath)
cdc$brkthruRatio = cdc$Vaccinated.with.outcome/(cdc$Vaccinated.with.outcome + cdc$Unvaccinated.with.outcome)

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

fl_rescale = 2148 # pop of FL/10k

# CURRENTLY ASSUMES SIM DATA STARTS ON WEDNESDAY, FEB 5, 2020
num_days       = nrow(d[d$serial==min(d$serial),])
num_replicates = length(unique(d$serial))
full_weeks     = floor((num_days-2)/7) # -2 because data starts on a Wednesday, want Friday
remainder      = (num_days-2) %% 7
week_idx       = c(c(0,0), rep(1:full_weeks, each=7), rep(full_weeks+1, remainder)) #
d$week_idx     = week_idx
weekly         = aggregate(d[,c('rcase', 'rdeath')], by=list(week_idx=d$week_idx, serial=d$serial), FUN=sum)

ed_tmp     = merge(ed, d[,c('date', 'week_idx')], by='date')
ed_weekly  = aggregate(ed_tmp[,c('rcase', 'rdeath')], by=list(week_idx=ed_tmp$week_idx), FUN=sum)
if (sum(ed_tmp$week_idx == min(ed_tmp$week_idx)) < 7) { ed_weekly = ed_weekly[-1,] }
if (sum(ed_tmp$week_idx == max(ed_tmp$week_idx)) < 7) { ed_weekly = ed_weekly[1:(nrow(ed_weekly)-1),] }
weekly$rdate = as.Date('2020-02-07') + (7 * 0:(length(unique(week_idx))-1))

names(ed_weekly) = c('week_idx', 'ercase', 'erdeath')
weekly = merge(weekly, ed_weekly, by='week_idx', all.x=T)

weekly_mean = aggregate(cbind(rcase,rdeath) ~ week_idx + rdate, data = weekly, mean, na.rm = TRUE)

d$wildtype = 1 - (d$vocprev1+d$vocprev2+d$vocprev3)

d_median = aggregate(cbind(wildtype,vocprev1,vocprev2,vocprev3,cinf,rcase,rdeath,inf,rhosp,VES,brkthruRatio,vaxInfs,unvaxInfs,hospInc,hospPrev,icuInc,icuPrev,vaxHosp,unvaxHosp,Rt,crcase,crdeath,VESAvg,brkthruRatioAvg) ~ date, data = d, mean, na.rm = TRUE)

panel_plot = function(filename, x, y, ymax, ylab, yaxt='s') {
    options(scipen=5)
    png(filename, width=1800, height=600, res=180)
    par(mar=c(3.1, 4.1, 0.5, 0.5), las=1)
    plot(x, y, type='n', xlab='', xaxt='n', yaxt=yaxt, ylab=ylab, xlim=c(min(x), as.Date('2022-05-01')), ylim=c(0,ymax), bty='n')
    title(xlab='Date (month)', line=2)
    shading()
    annotate('')
}

# Daily reported cases
ymax = max(d$rcase*fl_rescale, ed$rcase*fl_rescale, na.rm=T)
#ymax = 32894.47
panel_plot(paste0('projected_cases_220104', suffix, '.png'), d$date, d$rcase*fl_rescale, ymax, 'Reported cases per day', yaxt='n')
lines(ed$date, ed$rcase*fl_rescale, lwd=2)
for (serial in unique(d$serial)) {
    lines(d$date[d$serial == serial], d$rcase[d$serial == serial]*fl_rescale, col=adjustcolor('orange', alpha.f=0.04))
}
lines(d_median$date[d_median$date > '2021-12-01'], d_median$rcase[d_median$date > '2021-12-01']*fl_rescale, col=adjustcolor('orange', alpha.f=1), lwd=2)
yticks = seq(0,1e5,2e4)
axis(2, at=yticks, labels=paste0(yticks/1000, 'k'))
legend('topleft', legend=c('Reported cases', 'Modeled reported cases'), col=c('black', 'orange'), lty=1, lwd=3)
dev.off()

## Weekly reported cases
#ymax = max(weekly$rcase*fl_rescale, weekly$ercase*fl_rescale, na.rm=T)
#panel_plot('projected_cases-weekly_220104.png', weekly$rdate, weekly$rcase*fl_rescale, ymax, 'Reported cases per week')
#lines(weekly$rdate, weekly$ercase*fl_rescale, lwd=2)
#lines(weekly$rdate, weekly$rcase*fl_rescale, col='orange')
#lines(d_median$date[d_median$date > '2021-12-01'], d_median$rcase[d_median$date > '2021-12-01']*fl_rescale, col=adjustcolor('orange', alpha.f=1), lwd=2)
#legend('topleft', legend=c('Reported cases', 'Modeled reported cases'), col=c('black', 'orange'), lty=1, lwd=3)
#dev.off()

# Rt
panel_plot(paste0('projected_Rt_220104', suffix, '.png'), d$date, d$Rt, max(d$Rt), expression(R[t]))
#lines(d$date[1:(length(d$date)-9)], d$Rt[1:(length(d$date)-9)], col='red3')
for (serial in unique(d$serial)) {
    #lines(d$date[d$serial == serial], d$Rt[d$serial == serial], col=adjustcolor('red3', alpha.f=0.04))
    lines(d$date[d$serial == serial][1:(length(d$date)-9)], d$Rt[d$serial == serial][1:(length(d$date)-9)], col=adjustcolor('red3', alpha.f=0.04))
}
#lines(d_median$date[d_median$date > '2021-12-01'], d_median$Rt[d_median$date > '2021-12-01'], col=adjustcolor('red3', alpha.f=1), lwd=2)
lines(d_median$date[d_median$date > '2021-12-01' & d_median$date <= max(d_median$date) - 9], d_median$Rt[d_median$date > '2021-12-01' & d_median$date <= max(d_median$date) - 9], col=adjustcolor('red3', alpha.f=1), lwd=2)
abline(h=1.0, lty=3)
dev.off()

# All infections 
panel_plot(paste0('projected_infections_220104', suffix, '.png'), d$date, d$inf*fl_rescale, max(d$inf*fl_rescale), 'Modeled infections per day', yaxt='n')
#lines(d$date[-1], d$inf[-1]*fl_rescale, col='green4')
for (serial in unique(d$serial)) {
    lines(d$date[d$serial == serial], d$inf[d$serial == serial]*fl_rescale, col=adjustcolor('green4', alpha.f=0.04))
}
lines(d_median$date[d_median$date > '2021-12-01'], d_median$inf[d_median$date > '2021-12-01']*fl_rescale, col=adjustcolor('green4', alpha.f=1), lwd=2)
yticks = seq(0,1e6,2e5)
axis(2, at=yticks, labels=paste0(yticks/1000, 'k'))
dev.off()

# Reported deaths
ymax = max(d$rdeath*fl_rescale, ed$rdeath*fl_rescale, na.rm=T)
panel_plot(paste0('projected_deaths_220104', suffix, '.png'), d$date, d$rdeath*fl_rescale, ymax, 'Reported deaths per day')
lines(ed$date[1:(nrow(ed) - 14)], ed$rdeath[1:(nrow(ed) - 14)]*fl_rescale, lwd=2)
#lines(d$date, d$rdeath*fl_rescale, col='red')
for (serial in unique(d$serial)) {
    lines(d$date[d$serial == serial], d$rdeath[d$serial == serial]*fl_rescale, col=adjustcolor('red', alpha.f=0.02))
}
lines(d_median$date[d_median$date > '2021-12-01'], d_median$rdeath[d_median$date > '2021-12-01']*fl_rescale, col=adjustcolor('red', alpha.f=1), lwd=2)
legend('topleft', legend=c('Reported deaths', 'Modeled reported deaths'), col=c('black', 'red'), lty=1, lwd=3)
dev.off()

# Reported hospitalizations
#ymax = max(d$rhosp, ed$rhosp, na.rm=T)
#panel_plot('projected_hospitalizations_220104.png', d$date, d$rhosp, ymax, 'Reported hospitalizations')
#lines(ed$date, ed$rhosp, lwd=2)
#lines(d$date, d$rhosp, col='red')
#legend('topleft', legend=c('Reported hospitalizations', 'Modeled hospitalizations'), col=c('black', 'red'), lty=1, lwd=3)
#dev.off()

# VOC prevalence
panel_plot(paste0('projected_voc_prevalence_220104', suffix, '.png'), d$date, d$vocprev1, 1, 'Variant prevalence')
for (serial in unique(d$serial)) {
    lines(d$date[d$serial == serial], d$wildtype[d$serial == serial], col=adjustcolor('seagreen3', alpha.f=0.04))
    lines(d$date[d$serial == serial], d$vocprev1[d$serial == serial], col=adjustcolor('royalblue3', alpha.f=0.04))
    lines(d$date[d$serial == serial], d$vocprev2[d$serial == serial], col=adjustcolor('turquoise4', alpha.f=0.04))
    lines(d$date[d$serial == serial], d$vocprev3[d$serial == serial], col=adjustcolor('darkorchid3', alpha.f=0.04))
}
lines(d_median$date, d_median$wildtype, col='seagreen3')
lines(d_median$date, d_median$vocprev1, col='royalblue3')
lines(d_median$date, d_median$vocprev2, col='turquoise4')
lines(d_median$date, d_median$vocprev3, col='darkorchid3')
legend('topleft', legend=c('Wildtype', 'Alpha', 'Delta', 'Omicron'), col=c('seagreen3', 'royalblue3', 'turquoise4', 'darkorchid3'), lty=1, lwd=3)
dev.off()




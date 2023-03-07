#!/usr/bin/env Rscript
library(MMWRweek)

args = commandArgs(trailingOnly=TRUE)
# empty => -d .
# sniff for -d flag => process tail arg as a directory
# otherwise => process all args as files

if (length(args) == 0) {
  args = c("-d", ".")
}

if ("-d" %in% args) {
  srcfiles <- tail(args, 1) |> normalizePath() |> list.files(
    pattern = "plot_log[^\\]+\\.csv$", full.names = TRUE
  )
} else {
  srcfiles <- normalizePath(args)
}

covidabmloc <- getwd()
while(basename(covidabmloc) != "covid-abm") {
  if (dirname(covidabmloc) == covidabmloc) stop("not in covid-abm project: ", getwd())
  covidabmloc <- dirname(covidabmloc)
}
datapath <- file.path(covidabmloc, "exp", "active-vac")
refpath <- function(p) file.path(datapath, p)

calc_centered_avg = function(vals, halfwindow) {
    n       = length(vals)
    indices = 1:n
    ends    = pmin(indices + halfwindow + 1, n)
    starts  = ifelse(halfwindow >= indices, 0, indices - halfwindow)
    widths  = ends - starts
    ma      = mapply(function(s,e) sum(vals[s:e]), s = starts, e = ends) / widths
    return(ma)
}

#' @WARNING no logic here to select amongst these
#escambia_fraction    = 0.0153 # fraction of FL pop that lives in Escambia
pop_florida  = 21538187 #20609673
pop_escambia = 312212
pop_dade     = 2794464
per10k = 1e4/pop_florida

pop = read.table(refpath("../../pop/pseudo-300K/population-pseudo-300K.txt"), stringsAsFactors = F, header = T, sep = ' ')
bin_pops = as.integer(table(cut(pop$age, breaks = c(0, 5, 12, 18, 65, 120), right = F)))

ed = read.csv(refpath("rcasedeath-florida.csv"));
#ed = read.csv("rcasedeath-escambia.csv"); per10k = 1e4/pop_escambia
#ed = read.csv("rcasedeath-dade.csv"); per10k = 1e4/pop_dade
ed$date = as.Date(ed$Date)
#death_underreporting = 20.1/11.8 # excess death / recognized COVID deaths
ed$rcase = ed$rcase*per10k
#ed$rhosp = ed$rhosp*per10k
#ed$rdeath = ed$rdeath*per10k #*death_underreporting
#ed$rdeath = ed$death_incd*per10k
ed$rdeath  = ed$excess*per10k
ed$crcase  = cumsum(ed$rcase)
ed$crdeath = cumsum(ed$rdeath)

cdcFull = read.csv(refpath("Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Vaccination_Status.csv"), stringsAsFactors=F)
cdc = cdcFull[cdcFull$Age.group == "all_ages_adj" & cdcFull$Vaccine.product == "all_types"  & cdcFull$outcome == "case",]
cdc$MMWR.week_adj = as.numeric(sub('.{4}(.*)', '\\1', cdc$MMWR.week))
cdc$MMWR.year_adj = as.numeric(sub('(.{4}).*', '\\1', cdc$MMWR.week))
cdc$date = as.Date(MMWRweek2Date(cdc$MMWR.year_adj, cdc$MMWR.week_adj))
cdc$brkthruRatio = cdc$Vaccinated.with.outcome/(cdc$Vaccinated.with.outcome + cdc$Unvaccinated.with.outcome)
cdc$vaxOutcomeP10k = cdc$Vaccinated.with.outcome * (1e4/cdc$Fully.vaccinated.population)

hhsHosp = read.csv(refpath("COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries.csv"), stringsAsFactors=F, header=T)
hhsHosp = hhsHosp[hhsHosp$state == 'FL',]
hhsHosp$date = as.Date(hhsHosp$date)
hhsHosp = hhsHosp[order(hhsHosp$date),]
hhsHosp$hospInc = hhsHosp$previous_day_admission_adult_covid_confirmed*per10k

seroprev = read.csv(refpath("fl_seroprev_cdc.csv"), stringsAsFactors = F, header = T)
seroprev$date = as.Date(seroprev$date)

vax = read.table(refpath("state_based_counterfactual_doses.txt"), stringsAsFactors = F, header = T, sep = ' ')
vax = vax[vax$is_urg == 0 & vax$ref_location == "FL",]
vax$ndoses = 0
vax[vax$bin_min == 5,]$ndoses = vax[vax$bin_min == 5,]$n_doses_p10k * (bin_pops[2]/1e4)
vax[vax$bin_min == 12,]$ndoses = vax[vax$bin_min == 12,]$n_doses_p10k * (bin_pops[3]/1e4)
vax[vax$bin_min == 18,]$ndoses = vax[vax$bin_min == 18,]$n_doses_p10k * (bin_pops[4]/1e4)
vax[vax$bin_min == 65,]$ndoses = vax[vax$bin_min == 65,]$n_doses_p10k * (bin_pops[5]/1e4)
agg_vax = aggregate(x = list(doses = vax$ndoses), by = list(date=vax$date, dose=vax$dose), FUN = sum)
agg_vax$date = as.Date(agg_vax$date)
agg_vax$cov = agg_vax$doses/375474
cov_dates = unique(agg_vax$date)
cov1 = cumsum(agg_vax$cov[agg_vax$dose == 1])
cov2 = cumsum(agg_vax$cov[agg_vax$dose == 2])
cov3 = cumsum(agg_vax$cov[agg_vax$dose == 3])

doitall <- function(path) {
  fromdata = path; topng = file.path(
    dirname(path), basename(path) |>
      gsub(pattern = "plot_log([^\\]+)\\.csv$", replacement = "simvis\\1.png", x = _)
  )

  d = read.csv(fromdata)
  d$date = as.Date(d$date)

  d$crcase   = cumsum(d$rcase)
  d$crdeath  = cumsum(d$rdeath)
  is.na(d) = sapply(d, is.infinite)
  d$VESAvg = filter(d$VES, rep(1/7, 7), sides = 2)
  d$brkthruRatioAvg = filter(d$brkthruRatio, rep(1/7, 7), sides = 2)
  d$tot_std_doses = cumsum(d$std_doses)
  d$tot_urg_doses = cumsum(d$urg_doses)

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

  png(topng, height=1200, width=2400, res=200)
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

  # Seroprevalence
  ymax = max(d$cinf, d$seroprev, seroprev$upper, 1, na.rm=T)
  plot(d$date, d$cinf, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
  shading()
  abline(h=1, lty=3)
  arrows(seroprev$date, seroprev$all_seroprev_lower/100, seroprev$date, seroprev$all_seroprev_upper/100, length = 0.05, code=0, col = rgb(0,0,0,0.25))
  points(seroprev$date, seroprev$all_seroprev_point/100, pch = 20, cex = 0.5, col='darkorange')
  arrows(seroprev$date, seroprev$ped_seroprev_lower/100, seroprev$date, seroprev$ped_seroprev_upper/100, length = 0.05, code=0, col = rgb(0,0,0,0.25))
  points(seroprev$date, seroprev$ped_seroprev_point/100, pch = 20, cex = 0.5, col='royalblue')
  # lines(d$date, d$cinf, col='orangered', lty = 3)
  lines(d$date, d$seroprev, col='darkorange')
  lines(d$date, d$ped_seroprev, col='royalblue')
  annotate('Seroprevalence (total - orange; pediatric - blue)')

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
  ymax = max(cov1, d$cov1, 1, na.rm=T)
  plot(d$date, d$std_doses, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
  shading()
  abline(h=1, lty=3)
  lines(cov_dates, cov1, lty = 3)
  lines(cov_dates, cov2, lty = 2)
  lines(cov_dates, cov3, lty = 1)
  lines(d$date, d$cov1, col='purple', lty=3)
  lines(d$date, d$cov2, col='purple', lty=2)
  lines(d$date, d$cov3, col='purple', lty=1)
  annotate('Vaccination coverage (dose 1 - dotted; dose 2 - dashed; dose 3 = solid)')

  # hosp by vax status
  ymax = max(d$vaxHosp + d$unvaxHosp, rm.na=T)
  plot(d$date, d$vaxHosp, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
  shading()
  lines(hhsHosp$date, hhsHosp$hospInc)
  lines(d$date, d$vaxHosp, col='dodgerblue4')
  lines(d$date, d$unvaxHosp + d$vaxHosp, col='coral')
  annotate('Hosp by vax status (stacked)')

  # hosp inc/prev
  ymax = max(d$hospInc, d$hospPrev, rm.na=T)
  plot(d$date, d$hospPrev, type='n', xlab='', ylab='', xaxt='n', ylim=c(0,ymax), bty='n')
  shading()
  lines(d$date, d$hospPrev, col='orange')
  lines(hhsHosp$date, hhsHosp$hospInc)
  lines(d$date, d$hospInc, col='orange', lty=3)
  annotate('Hosp prev (solid) and inc (dotted)')

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
}

sapply(srcfiles, doitall)

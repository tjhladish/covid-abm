setwd("~/documents/work/covid-abm/exp/vac_impact")

library(dplyr)

# READING IN EMPIRICAL FL VACCINATION DATA FROM CDC
fl_dat = read.table(file = "../../pop/fl_vac/trends_in_number_of_covid19_vaccinations_in_fl.csv", header = T, stringsAsFactors = F, sep = ',', skip = 2)
fl_cut = fl_dat %>% 
  filter(Date.Type == "Admin" & Date != "N/A") %>%
  select(date_str = Date, daily_dose1_plus = Daily.Count.People.Receiving.Dose.1, cum_dose1_plus = People.Receiving.1.or.More.Doses.Cumulative)
fl_cut$date = as.numeric(as.Date(fl_cut$date_str)-as.Date("2020-02-05"))
fl_cut$cum_sum_daily = cumsum(fl_cut$daily_dose1_plus)

# FL DOH VAX DATA
vax_dat = read.table(file = "../../pop/fl_vac/fl_vac_v3.txt", header = T, stringsAsFactors = F, sep = ' ')
vax_adj = vax_dat %>% group_by(date) %>% summarise(dose_sum = sum(mrna_first_doses, mrna_sec_doses, jj_doses)) %>% ungroup()

# DATA FROM BEFORE BUGFIX
raw_dat = read.table(file = "./vac_sch_debug.csv", header = F, stringsAsFactors = F, sep = ',')
dat = data.frame(date = raw_dat$V6, bin_min = raw_dat$V7, n = raw_dat$V8, p = raw_dat$V9, num_sch = raw_dat$V10)

daily_sum_sch = dat %>% group_by(date) %>% arrange(date) %>% summarise(sum_sch = sum(num_sch)) %>% ungroup()
daily_sum_sch$cum_sch = cumsum(daily_sum_sch$sum_sch)
plot(x = daily_sum_sch$date, y = daily_sum_sch$sum_sch, type = "l")
plot(x = daily_sum_sch$date, y = daily_sum_sch$cum_sch, type = "l", col = "red")

# DATA FROM AFTER BUGFIX
raw_dat2 = read.table(file = "./vac_sch_debug2.csv", header = F, stringsAsFactors = F, sep = ',')
dat2 = data.frame(date = raw_dat2$V6, bin_min = raw_dat2$V7, n = raw_dat2$V8, p = raw_dat2$V9, num_sch = raw_dat2$V10)

daily_sum_sch2 = dat2 %>% group_by(date) %>% arrange(date) %>% summarise(sum_sch = sum(num_sch)) %>% ungroup()
daily_sum_sch2$cum_sch = cumsum(daily_sum_sch2$sum_sch)
plot(x = daily_sum_sch2$date, y = daily_sum_sch2$sum_sch, type = "l")
plot(x = daily_sum_sch2$date, y = daily_sum_sch2$cum_sch, type = "l", col = "red")

fl_vax_dose1plus = vax_dat %>% group_by(date) %>% summarise(dose1plus = sum(mrna_first_doses, jj_doses)) %>% ungroup()
fl_vax_dose1plus$cum_dose1plus = cumsum(fl_vax_dose1plus$dose1plus)

png(filename = "vaccine_uptake.png", width = 1600, height = 1200, res = 180)
plot(x = as.Date("2020-02-05")+fl_cut$date, y = fl_cut$cum_sum_daily/21538187, type = "l", xlab = "Date", ylab = "Vaccination Coverage (1+ Dose)") # FL cum dose1+ data (adj. for pop ratio)
#lines(x = as.Date("2020-02-05")+day_sum_sch_dose$date, y = day_sum_sch_dose$cum_non_hcw, type = "l", col = "orange")         # model cum dose1 data
lines(x = as.Date("2020-02-05")+day_sum_sch_dose$date, y = day_sum_sch_dose$cum_sch/375474, type = "l", col = "red")         # model cum dose1 data
#lines(x = as.Date("2020-02-05")+day_sum_sch_dose$date, y = day_sum_sch_dose$cum_doses, type = "l", col = "green")         # model cum dose1 data
#lines(x = as.Date("2020-02-05")+daily_sum_sch2$date, y = daily_sum_sch2$cum_sch+29185, type = "l", col = "green") # model adj. for vaccinations before model data
points(x = as.Date(fl_vax_dose1plus$date), y = fl_vax_dose1plus$cum_dose1plus/21538187)
abline(h = 0, lty = 3)
legend('topleft', legend = c("CDC Data", "FDOH Data", "Model"), col = c("black", "black", "red"), lty = c(1, NA, 1), pch = c(NA, 1, NA), bty = "n")
dev.off()

# CHECKING DOSES ALLOCATED VS. PPL SCHEDULED
daily_sch_dose = read.table(file = "./daily_sch_dose.out", header = F, sep = ' ', stringsAsFactors = F)
colnames(daily_sch_dose) = c("date", "sch", "non_hcw", "doses")
day_sum_sch_dose = daily_sch_dose %>% group_by(date) %>% summarise(sum_sch = sum(sch), sum_doses = sum(doses), max_doses = max(doses), sum_non_hcw = sum(non_hcw)) %>% ungroup()
day_sum_sch_dose$cum_sch = cumsum(day_sum_sch_dose$sum_sch)
day_sum_sch_dose$cum_non_hcw = cumsum(day_sum_sch_dose$sum_non_hcw)
day_sum_sch_dose$cum_doses = cumsum(day_sum_sch_dose$max_doses)
plot(x = as.Date("2020-02-05")+day_sum_sch_dose$date, y = day_sum_sch_dose$max_doses, type = "l", col = "red", ylim = c(0, max(day_sum_sch_dose$max_doses)))     # doses available
lines(x = as.Date("2020-02-05")+day_sum_sch_dose$date, y = day_sum_sch_dose$sum_doses/2, type = "l", col = "blue") # 1/2 doses available
lines(x = as.Date("2020-02-05")+day_sum_sch_dose$date, y = day_sum_sch_dose$sum_sch, type = "l")                   # cum ppl sch
lines(x = as.Date(vax_adj$date), y = vax_adj$dose_sum*(375474/21538187)/7, type = "l", col = "green")


png(filename = "doses_per_day.png", width = 1600, height = 1200, res = 180)
plot(x = as.Date("2020-02-05")+day_sum_sch_dose$date, y = day_sum_sch_dose$max_doses/37.5474, type = "l", col = "red", ylim = c(0,150),
     xlab = "Date", ylab = "Doses Administered per 10K")     # doses available
points(x = as.Date(vax_adj$date), y = vax_adj$dose_sum*(10000/21538187)/7)
legend('topright', legend = c("FDOH Data", "Model"), col = c("black", "red"), lty = c(NA, 1), pch = c(1, NA), bty = "n")
dev.off()
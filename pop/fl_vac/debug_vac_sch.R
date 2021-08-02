setwd("~/documents/work/covid-abm/exp/vac_impact")
library(dplyr)

### USEFUL VARIABLES
sim_start_date = as.Date("2020-02-05") # simulation start date
fl_pop = 21538187                      # population of Florida (2020 census)
sim_pop = 375474                       # 300K synthetic population size

### READING IN EMPIRICAL FL VACCINATION DATA FROM CDC
fl_dat = read.table(file = "../../pop/fl_vac/trends_in_number_of_covid19_vaccinations_in_fl.csv", header = T, stringsAsFactors = F, sep = ',', skip = 2)
# filter out unwanted/empty rows; select and rename desired fields
fl_cut = fl_dat %>% filter(Date.Type == "Admin" & Date != "N/A") %>%
  select(date = Date, daily_dose1_plus = Daily.Count.People.Receiving.Dose.1, cum_dose1_plus = People.Receiving.1.or.More.Doses.Cumulative)
# calculate simulation day from date
fl_cut$sim_day = as.numeric(as.Date(fl_cut$date)-sim_start_date)
# cast date as proper type
fl_cut$date = as.Date(fl_cut$date)
# calculate cumulative sum of daily doses (1+)
fl_cut$cum_sum_daily = cumsum(fl_cut$daily_dose1_plus)

### FL DOH VAX DATA
vax_dat = read.table(file = "../../pop/fl_vac/fl_vac_v3.txt", header = T, stringsAsFactors = F, sep = ' ')
# calculate sum of doses administered per week
vax_adj = vax_dat %>% group_by(date) %>% summarise(dose_sum = sum(mrna_first_doses, mrna_sec_doses, jj_doses)) %>% ungroup()
# cast date as proper type
vax_adj$date = as.Date(vax_adj$date)

### DATA FROM BEFORE BUGFIX
raw_dat = read.table(file = "./vac_sch_debug.csv", header = F, stringsAsFactors = F, sep = ',')
# select and rename desired fields
dat = data.frame(sim_day = raw_dat$V6, bin_min = raw_dat$V7, n = raw_dat$V8, p = raw_dat$V9, num_sch = raw_dat$V10)
# calculate sum of people scheduled each simulated day
daily_sum_sch = dat %>% group_by(sim_day) %>% arrange(sim_day) %>% summarise(sum_sch = sum(num_sch)) %>% ungroup()
# calculate cumulative sum of people scheduled
daily_sum_sch$cum_sch = cumsum(daily_sum_sch$sum_sch)

### DATA FROM AFTER BUGFIX
raw_dat2 = read.table(file = "./vac_sch_debug2.csv", header = F, stringsAsFactors = F, sep = ',')
# select and rename desired fields
dat2 = data.frame(sim_day = raw_dat2$V6, bin_min = raw_dat2$V7, n = raw_dat2$V8, p = raw_dat2$V9, num_sch = raw_dat2$V10)
# calculate sum of people scheduled each simulated day
daily_sum_sch2 = dat2 %>% group_by(sim_day) %>% arrange(sim_day) %>% summarise(sum_sch = sum(num_sch)) %>% ungroup()
# calculate cumulative sum of people scheduled
daily_sum_sch2$cum_sch = cumsum(daily_sum_sch2$sum_sch)

### CHECKING DOSES ALLOCATED VS. PPL SCHEDULED
daily_sch_dose = read.table(file = "./daily_sch_dose.out", header = F, sep = ' ', stringsAsFactors = F)
# rename columns
colnames(daily_sch_dose) = c("sim_day", "sch", "non_hcw", "doses")
# calculate desired variables for each simulated day: sum of people (all and non-hcw) scheduled and maximum doses available
day_sum_sch_dose = daily_sch_dose %>% group_by(sim_day) %>%
  summarise(sum_sch = sum(sch), max_doses = max(doses), sum_non_hcw = sum(non_hcw)) %>% ungroup()
# calculate date from sim_day
day_sum_sch_dose$date = sim_start_date+day_sum_sch_dose$sim_day
# calculate cumulative sum of scheduled people
day_sum_sch_dose$cum_sch = cumsum(day_sum_sch_dose$sum_sch)
# calculate cumulative sum of non-hcw people scheduled
day_sum_sch_dose$cum_non_hcw = cumsum(day_sum_sch_dose$sum_non_hcw)
# calculate cumulative sum of daily doses available
day_sum_sch_dose$cum_doses = cumsum(day_sum_sch_dose$max_doses)

### FLDOH DATA ADDITIONAL PROCESSING
fl_vax_dose1plus = vax_dat %>% group_by(date) %>% summarise(dose1plus = sum(mrna_first_doses, jj_doses)) %>% ungroup()
# cast date as proper type
fl_vax_dose1plus$date = as.Date(fl_vax_dose1plus$date)
fl_vax_dose1plus$cum_dose1plus = cumsum(fl_vax_dose1plus$dose1plus)

### PLOTTING VARIABLES
png_width = 1600
png_height = 1200
png_res = 180

### PLOTTING
# pre-bugfix plots of daily and cumulative people scheduled
plot(x = daily_sum_sch$sim_day, y = daily_sum_sch$sum_sch, type = "l")
plot(x = daily_sum_sch$sim_day, y = daily_sum_sch$cum_sch, type = "l", col = "red")

# post-bugfix plots of daily and cumulative people scheduled
plot(x = daily_sum_sch2$sim_day, y = daily_sum_sch2$sum_sch, type = "l")
plot(x = daily_sum_sch2$sim_day, y = daily_sum_sch2$cum_sch, type = "l", col = "red")

# plot comparing CDC, FLDOH and model data on cumulative vaccination coverage (1+ dose)
png(filename = "vaccine_uptake.png", width = png_width, height = png_height, res = png_res)
plot(x = fl_cut$date, y = fl_cut$cum_sum_daily/fl_pop, type = "l", xlab = "Date", ylab = "Vaccination Coverage (1+ Dose)") # FL cum dose1+ data (adj. for pop ratio)
lines(x = day_sum_sch_dose$date, y = day_sum_sch_dose$cum_sch/sim_pop, type = "l", col = "red")                            # model cum dose1 data
points(x = fl_vax_dose1plus$date, y = fl_vax_dose1plus$cum_dose1plus/fl_pop)
abline(h = 0, lty = 3)
legend('topleft', legend = c("CDC Data", "FDOH Data", "Model"), col = c("black", "black", "red"), lty = c(1, NA, 1), pch = c(NA, 1, NA), bty = "n")
dev.off()

# plot comparing FLDOH and model data on doses administered per 10K
png(filename = "doses_per_day.png", width = png_width, height = png_height, res = png_res)
plot(x = day_sum_sch_dose$date, y = day_sum_sch_dose$max_doses/(sim_pop/10000), type = "l", col = "red", ylim = c(0,150),
     xlab = "Date", ylab = "Doses Administered per 10K")     # doses available
points(x = vax_adj$date, y = vax_adj$dose_sum*(10000/fl_pop)/7)
legend('topright', legend = c("FDOH Data", "Model"), col = c("black", "red"), lty = c(NA, 1), pch = c(1, NA), bty = "n")
dev.off()
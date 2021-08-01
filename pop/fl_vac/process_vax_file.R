setwd("~/documents/work/covid-abm/pop/fl_vac/")

library(dplyr)

dat = read.table(file = "./fl-vac-dat_20210123-0529.csv", header =  T, stringsAsFactors = F, sep = ",")
dat$date = as.Date(format(dat$date, format = "%Y%m%d"))

fl_dat = dat %>% filter(County == "Florida") %>% group_by(date, age_grp) %>% arrange(date, desc(age_grp)) %>% ungroup()
fl_dat$age_grp = gsub(" ", "_", fl_dat$age_grp)

tmp_bin_dat = list()
for(bin in unique(fl_dat$age_grp)) {
  tmp_bin_dat[[bin]] = fl_dat %>% filter(age_grp == bin)
}

fl_dat_diff = data.frame(age_grp = character(), first_dose = integer(), complete = integer(), vac_pers = integer(), date = as.Date(character()))
for(bin in unique(fl_dat$age_grp)) {
  print(bin)
  for(row in 1:nrow(tmp_bin_dat[[bin]])) {
    tmp_age = tmp_bin_dat[[bin]]$age_grp[row]
    tmp_date = tmp_bin_dat[[bin]]$date[row]
    if(tmp_bin_dat[[bin]]$date[row] == min(tmp_bin_dat[[bin]]$date)) {
      fd_diff = tmp_bin_dat[[bin]]$first_dose[row]
      c_diff = tmp_bin_dat[[bin]]$complete[row]
      vp_diff = tmp_bin_dat[[bin]]$vac_pers[row]
    } 
    else if(tmp_bin_dat[[bin]]$date[row] == tmp_bin_dat[[bin]]$date[row-1]+7) {
      fd_diff = tmp_bin_dat[[bin]]$first_dose[row]-tmp_bin_dat[[bin]]$first_dose[row-1]
      c_diff = tmp_bin_dat[[bin]]$complete[row]-tmp_bin_dat[[bin]]$complete[row-1]
      vp_diff = tmp_bin_dat[[bin]]$vac_pers[row]-tmp_bin_dat[[bin]]$vac_pers[row-1]
    }
    
    tmp_df = data.frame(age_grp = tmp_age, first_dose = fd_diff, complete = c_diff, vac_pers = vp_diff, date = tmp_date)
    fl_dat_diff = fl_dat_diff %>% add_row(tmp_df)
  }
}

JJprop = 0.08027
fl_vac_adjusted = data.frame(age_grp = fl_dat_diff$age_grp, date = fl_dat_diff$date,
                             mrna_first_doses = (fl_dat_diff$vac_pers - (round(JJprop * fl_dat_diff$complete))), 
                             mrna_sec_doses = (fl_dat_diff$complete - (round(JJprop * fl_dat_diff$complete))), 
                             jj_doses = round(JJprop * fl_dat_diff$complete))

final_dat = fl_vac_adjusted %>% ungroup() %>% group_by(date, age_grp) %>% arrange(date, age_grp) %>% ungroup()

final_dat$bin_min = as.numeric(substring(final_dat$age_grp, 1,2))

final_dat$bin_max = substring(final_dat$age_grp, 4,5)
final_dat$bin_max[final_dat$bin_max == "_y"] = "120"
final_dat$bin_max = as.numeric(final_dat$bin_max)

write.table(final_dat, file = "fl_vac_v3.txt", row.names = F, sep = ' ', quote = F)

### LOOKING AT CUMULATIVE DOSE SUMS
dose_sum = final_dat %>% group_by(date) %>% summarise(mrna_1_sum = sum(mrna_first_doses), mrna_2_sum = sum(mrna_sec_doses), jj_sum = sum(jj_doses)) %>% ungroup()

dose_sum$cum_mrna_1_sum = cumsum(dose_sum$mrna_1_sum)
dose_sum$cum_mrna_2_sum = cumsum(dose_sum$mrna_2_sum)
dose_sum$cum_jj_sum = cumsum(dose_sum$jj_sum)
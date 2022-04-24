library(DBI)
library(data.table)
library(lubridate)
library(ggplot2)
library(viridis)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/fl_vac_counterfactual/") }

.args <- if (interactive()) c(
  "sim_data_0.sqlite"
) else commandArgs(trailingOnly = TRUE)

db_path <- .args[1]

serial <- sub("sim_data_([[:digit:]]*).sqlite", "\\1", db_path)
figPath <- paste0("./fig/sim_", serial, "_diagnostics")
dir.create(figPath, showWarnings = FALSE)

# OFFSPRING DISTRIBUTION
db = dbConnect(RSQLite::SQLite(), db_path)
stmt = "select * from infection_history;"
seconday_inf_counts = dbGetQuery(db, stmt)
dbDisconnect(db)

seconday_inf_counts = setDT(seconday_inf_counts)

hist(seconday_inf_counts[, num_sec_infs], freq=F)

out_table = seconday_inf_counts[, .N, by = .(num_sec_infs)][order(num_sec_infs)]
fwrite(x=out_table, file='./offspring_distribution.csv', sep=',', quote=F)

out_table = seconday_inf_counts[, .N, by = .(day=infected_time, num_sec_infs)][order(day, num_sec_infs)]
fwrite(x=out_table, file='./daily_offspring_distributions.csv', sep=',', quote=F)

# INFECTIONS BY TYPE
db = dbConnect(RSQLite::SQLite(), db_path)
db2 = dbConnect(RSQLite::SQLite(), "../../pop/sim_pop-pseudo-300K-3.1/sim_pop-pseudo-300K-3.1.sqlite")
dbExecute(db, "attach database '../../pop/sim_pop-pseudo-300K-3.1/sim_pop-pseudo-300K-3.1.sqlite' as synthpop")

stmt = "select ih.inf, ih.infected_time, ih.inf_place_id, l.type as inf_place_type, ih.inf_owner_id, p.age, r.locid as home_id, m.locid as day_id, m.type as day_type 
from infection_history as ih 
left join synthpop.reside as r on (ih.inf_owner_id + 1) = r.pid 
left join synthpop.movement as m on (ih.inf_owner_id + 1) = m.pid
left join synthpop.loc as l on (ih.inf_place_id + 1) = l.locid
left join synthpop.pers as p on (ih.inf_owner_id + 1) = p.pid;"

inf_hist_w_locs = dbGetQuery(db, stmt)
dbDisconnect(db)
dbDisconnect(db2)

inf_hist_w_locs = setDT(inf_hist_w_locs)
inf_hist_w_locs[, inf_place_id_adj := inf_place_id + 1]
inf_hist_w_locs[is.na(day_id), day_id := -1]

assign_specifier <- function(SD) {
  specifier = ''
  if (SD[, inf_place_id_adj] == 0) { specifier = ''}
  switch(SD[, inf_place_type],
         w = {
           if (SD[day_id == inf_place_id_adj, .N]) { specifier = 'work_staff'}
           else { specifier = 'work_patron' }
         },
         h = {
           if (SD[home_id == inf_place_id_adj, .N]) { specifier = 'in_home'}
           else { specifier = 'social' }
         },
         s = {
           if (SD[age >= 18, .N]) { specifier = 'school_staff'}
           else { specifier = 'student' }
         },
         hf = {
           if (SD[day_id == inf_place_id_adj, .N]) { specifier = 'hcw'}
           else { specifier = 'patient' }
         },
         n = {
           if (SD[day_id == inf_place_id_adj, .N]) { specifier = 'ltcf_staff'}
           else { specifier = 'ltcf_resident' }
         })
  return(specifier)
}

inf_hist_w_locs[, inf_specifier := assign_specifier(.SD), by = row.names(inf_hist_w_locs)]

inf_by_loc_ts = inf_hist_w_locs[inf_place_id != -1, .N, by = .(infected_time, inf_place_type, inf_specifier)][order(infected_time)]

sim_start_date = ymd("2020-02=05")
loc_labs = c("Home", "Work", "Hospital", "School", "LTCF")
names(loc_labs) = c("h", "w", "hf", "s", "n")

inf_by_loc_plot = ggplot(inf_by_loc_ts, aes(x = sim_start_date + infected_time, y = N)) +
  geom_line(aes(group = inf_specifier, color = inf_specifier)) +
  facet_grid(rows = vars(inf_place_type), scales = "free", labeller = labeller(inf_place_type = loc_labs)) +
  scale_color_viridis_d(name = "", 
                        limits = c("in_home", "social", "work_staff", "work_patron", "school_staff", "student", "ltcf_staff", "hcw", "ltcf_resident", "patient"),
                        labels = c("In-household", "Social", "Work staff", "Business patron", "School staff", "Student", "LTCF staff", "HCW", "LTCF resident", "Patient")) +
  xlab("Date") +
  ylab("Num. infections")

ggsave(filename = paste0(figPath, "/inf_by_loc.png"), plot = inf_by_loc_plot, device = 'png', units = 'in', height = 6, width = 12, dpi = 300)
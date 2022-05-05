library(DBI)
library(data.table)
library(lubridate)
library(ggplot2)
library(viridis)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/flvtms-vac_comp-300k") }

.args <- if (interactive()) c(
  "sim_data_100.sqlite",
  "../../pop/sim_pop-pseudo-300K-3.1/sim_pop-pseudo-300K-3.1.sqlite"
) else commandArgs(trailingOnly = TRUE)

db_path <- .args[1]

serial <- sub("sim_data_([[:digit:]]*).sqlite", "\\1", db_path)
figPath <- paste0("./fig/sim_", serial, "_diagnostics")
dir.create(figPath, showWarnings = FALSE)

# OFFSPRING DISTRIBUTION
db = dbConnect(RSQLite::SQLite(), db_path)
stmt = "select * from infection_history;"
inf_history = dbGetQuery(db, stmt)
dbDisconnect(db)

inf_history = setDT(inf_history)

hist(inf_history[, num_sec_infs], freq=F)

out_table = inf_history[, .N, by = .(num_sec_infs)][order(num_sec_infs)]
ofilename = paste0("./offspring_distribution_", serial, ".csv")
fwrite(x=out_table, file=ofilename, sep=',', quote=F)

out_table = inf_history[, .N, by = .(day=infected_time, num_sec_infs)][order(day, num_sec_infs)]
ofilename = paste0("./daily_offspring_distributions_", serial, ".csv")
fwrite(x=out_table, file=ofilename, sep=',', quote=F)

# INFECTIONS BY TYPE
db = dbConnect(RSQLite::SQLite(), db_path)
stmt = paste0("attach database '", .args[2], "' as synthpop")
dbExecute(db, stmt)

stmt = "select ih.inf, ih.infected_time, ih.inf_place_id, l.type as inf_place_type, ih.inf_owner_id, p.age, r.locid as home_id, m.locid as day_id, m.type as day_type 
from infection_history as ih 
left join synthpop.reside as r on (ih.inf_owner_id + 1) = r.pid 
left join synthpop.movement as m on (ih.inf_owner_id + 1) = m.pid
left join synthpop.loc as l on (ih.inf_place_id + 1) = l.locid
left join synthpop.pers as p on (ih.inf_owner_id + 1) = p.pid;"

inf_hist_w_locs = dbGetQuery(db, stmt)
dbDisconnect(db)

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

# AVG TIMING BETWEEN INFECTION AND DEATH
db = dbConnect(RSQLite::SQLite(), db_path)
stmt = 'select ih.inf, ih.inf_owner_id, id.detected_state, id.reported_time, ih.death_time, (ih.death_time - id.reported_time) as lag 
from infection_history as ih  
left join infection_detection as id on (ih.inf = id.inf) 
where ih.detected != 0 and ih.death_time != 2147483647;'
det_to_dth_lag = dbGetQuery(db, stmt)
det_to_dth_lag = setDT(det_to_dth_lag)
dbDisconnect(db)

hist(det_to_dth_lag$lag, freq = F)
hist(det_to_dth_lag[lag >= 0, lag], freq = F)

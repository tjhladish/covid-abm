library(ggplot2)
library(data.table)
library(patchwork)
library(lubridate)
library(ggsci)
library(DBI)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/active-vac/") }

.args <- if (interactive()) c(
  "covid-active-v6.1.sqlite",
  "v6-vaclogs"
) else commandArgs(trailingOnly = TRUE)

if (length(.args) != 2) {
  stop("Usage: [path to batch database] [path to dir with vpdr files (no trailing /)]")
}

get_alloc = function(alloc_code) {
  alloc = switch (as.character(alloc_code),
    "1" = "LIC",
    "2" = "MIC",
    "3" = "HIC",
    "4" = "USA"
  )
  return(alloc)
}

get_scenario = function(db_path, sim_serial) {
  select_query = paste0("select * from par where serial = ", sim_serial, ";")
  con = dbConnect(RSQLite::SQLite(), db_path)
  scen_query = dbGetQuery(con, select_query) |> setDT()
  dbDisconnect(con)
  
  if (scen_query[, pas_vac] == 1) {
    strat = "Passive"
    alloc = get_alloc(scen_query[, pas_alloc])
  } else {
    strat = switch (as.character(scen_query[, act_vac]),
                    "1" = "Ring",
                    "2" = "Risk (hosp.)",
                    "3" = "Risk (age)"
    )
    alloc = get_alloc(scen_query[, act_alloc])
  }
  
  quar = ifelse((scen_query[, quar] == 0), "no_quar", "quar")
  inf_cond = ifelse((scen_query[, inf_con] == 2), "non_case_only", "any_status")
  
  return(c("vax_strat" = strat, "vax_alloc" = alloc, "quar" = quar, "inf_cond" = inf_cond))
}

age_vaxd_by_day = function(db_path, file_path) {
  serial = sub("v6-vaclogs/vpdr([[:digit:]]*).csv", "\\1", file_path)
  scenario = get_scenario(db_path, serial)
  
  d = fread(file_path)
  d[avg_age_dose1 == -1, avg_age_dose1 := NA]
  d[avg_age_dose2 == -1, avg_age_dose2 := NA]
  d[avg_age_dose3 == -1, avg_age_dose3 := NA]
  
  out_d = melt(d[, .(serial, day, "1" = avg_age_dose1, "2" = avg_age_dose2, "3" = avg_age_dose3)],
               id.vars = c("serial", "day"),
               variable.name = "dose",
               value.name = "avg_age")
  out_d[, `:=`(
    vax_strat = as.character(scenario["vax_strat"]),
    vax_alloc = as.character(scenario["vax_alloc"]),
    quar = as.character(scenario["quar"]),
    inf_cond = as.character(scenario["inf_cond"])
  )]
  
  
  return(out_d)
}

#' sims to run
#' 724000 --> ring (1), USA, no quar, no inf constraint
#' 728000 --> risk hosp (2), USA, no quar, no inf constraint
#' 732000 --> risk age (3), USA, no quar, no inf constraint
#' 466000 --> passive, USA, no quar, no inf constraint

final_dt = rbindlist(lapply(list.files(.args[2], full.names = TRUE), age_vaxd_by_day, db_path = .args[1]))

p = ggplot(final_dt) +
  geom_line(aes(x = day, y = avg_age, color = vax_strat), linewidth = 1, alpha = 0.5) +
  facet_grid(cols = vars(paste0("Dose ", dose))) +
  scale_color_aaas(name = element_blank()) +
  labs(x = element_blank(), y = "Median age") +
  theme_light() +
  theme(legend.position = "bottom", text = element_text(size = 15))

ggsave("avg_age_vaxd.png",
  plot = p,
  device = png,
  width = 10,
  height = 5,
  units = "in",
  dpi = 200
)
library(ggplot2)
library(data.table)
library(patchwork)
library(lubridate)
library(ggsci)
library(DBI)
library(scales)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/active-vac/") }
load("./fig/vis_support.rda")

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
    strat = "Standard"
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
final_dt[, seven_rolling_mean_avg_age := filter(avg_age, rep(1/7, 7), sides = 2), by = .(serial)]

iqr_dt = final_dt[, .(
  median = as.numeric(quantile(avg_age, na.rm = TRUE, probs = c(0.5))),
  upper = as.numeric(quantile(avg_age, na.rm = TRUE, probs = c(0.95))),
  lower = as.numeric(quantile(avg_age, na.rm = TRUE, probs = c(0.05)))
), by = .(day, dose, vax_strat, vax_alloc, quar, inf_cond)]

sim_start_date = ymd("2020-02-10")
iqr_dt[, date := sim_start_date + day]

text_labs = data.table(
  x = ymd("2021-01-05"),
  y = c(5, 44),
  lab = c("Minimum vaccinatable age", "Average vaccinatable population age"),
  dose = factor(3, levels = c("1","2","3"))
)

p = ggplot(iqr_dt[date >= ymd("2021-01-01") & date <= ymd("2022-03-31")]) +
  geom_month_background(
    iqr_dt[date >= ymd("2021-01-01") & date <= ymd("2022-03-31")], 
    by = c("dose"),
    value = "upper", 
    ymin = 0, 
    ymax = 100,
    font.size = 5
  ) +
  geom_ribbon(aes(x = date, ymin = lower, ymax = upper, fill = vax_strat), alpha = 0.25) +
  geom_line(aes(x = date, y = median, color = vax_strat), linewidth = 1, alpha = 1) +
  geom_hline(aes(yintercept = 5), linetype = "42") +
  geom_hline(aes(yintercept = 44), linetype = "42") +
  geom_text(
    data = text_labs,
    aes(x = x, y = y, label = lab),
    size = 3,
    hjust = "bottom",
    vjust = "left",
    nudge_y = 1
  ) +
  facet_grid(cols = vars(paste0("Dose ", dose))) +
  scale_color_aaas(name = element_blank()) +
  scale_fill_aaas(name = element_blank()) +
  labs(x = element_blank(), y = "Mean age") +
  theme_minimal() +
  scale_x_null() +
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
  coord_cartesian(expand = FALSE, clip = "off") +
  # ylim(c(0, 100)) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(2, "lines")
  )

ggsave("avg_age_vaxd.png",
  plot = p,
  device = png,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

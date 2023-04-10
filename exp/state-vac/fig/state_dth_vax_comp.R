library(data.table)
library(ggplot2)
library(scales)

#' ACS_2019_pop_data.csv available in active vac dir
#' cdc_covid-19_vax_data.csv available from the download_cdc_covid_vax_data.sh script
#' United_States_COVID-19_Cases_and_Deaths_by_State_over_Time.csv available from https://data.cdc.gov/Case-Surveillance/United-States-COVID-19-Cases-and-Deaths-by-State-o/9mfq-cb36
#' COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries.csv available from https://healthdata.gov/Hospital/COVID-19-Reported-Patient-Impact-and-Hospital-Capa/g62h-syeh/data
.args <- if (interactive()) c(
  file.path("../../raw_data/process/cdc_vax.rds"),
  file.path("../../raw_data/process/us_outcomes.rds"),
  file.path("fig/vis_support.rda"),
  file.path("fig/output/state_dth_vax_comp.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[3])

vax.dt <- readRDS(.args[1])[dose == 1][
  order(date), .(
    date, outcome = "dose1",
    value = cumsum(tot_prop)
  ), by=.(state = location)
]

out.dt <- readRDS(.args[2])[, .(state, date, outcome, value = cinc_p10k)]

plt.dt <- rbind(vax.dt, out.dt)

plt.dt[, outcome := factor(outcome, levels = c("dose1", "hosp", "death"), ordered = TRUE)]

p <- ggplot(plt.dt[between(date, "2021-03-01", "2022-08-31")]) + aes(date, value, color = state) +
  facet_grid(
    outcome ~ ., scales = "free_y", switch = "y",
    labeller = labeller(outcome = c(
      dose1 = "Dose 1 Coverage,\n% of Population",
      hosp = "Cumulative Hospitalizations,\nPer 10k Population",
      death = "Cumulative Deaths,\nPer 10k Population"
    ))
  ) +
  geom_month_background(
    plt.dt[between(date, "2021-03-01", "2022-08-31")],
    by = c("outcome"),
    font.size = 3
  ) +
  geom_line() +
  coord_cartesian(expand = FALSE, clip = "off") +
  scale_color_discrete(NULL) +
  theme_minimal() +
  theme(
    strip.placement = "outside",
    axis.title = element_blank(),
    panel.spacing.y = unit(1, "line"),
    legend.position = "bottom",
    axis.text.x = element_blank()
  )

ggsave(tail(.args, 1), p, units = "px", width = 3600, height = 3600, dpi = 360)

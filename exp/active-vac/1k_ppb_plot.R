library(data.table)
library(lubridate)
library(ggplot2)
library(dplyr)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/active-vac/") }

.args <- if (interactive()) c(
  "ppb_fits.csv"
) else commandArgs(trailingOnly = TRUE)

output_path = file.path('.', 'fig')
dir.create(output_path, recursive = TRUE)

d = fread(.args[1])
melt.d = melt(d, id.vars = c("date"), variable.name = "serial", value.name = "anchor_val")

mean.d = melt.d[, .(mean_val = mean(anchor_val)), by = .(date)]
fwrite(x = mean.d, file = "./1k_mean_ppb.csv", sep = ",")

label_first_mth_and_januarys <- function(x) {
  ifelse(is.na(lag(x)) | year(lag(x)) != year(x),
         paste0(month(x, label = T), "\n", year(x)),
         paste0(month(x, label = T))
  )
}

ppb_plot = ggplot(melt.d) +
  geom_line(aes(x = date, y = anchor_val, group = serial), alpha = 0.01) +
  theme_light() +
  theme(legend.position = 'none', text = element_text(size = 20)) +
  ylim(c(0,1)) +
  labs(x = "Date", y = "PPB value") +
  scale_x_date(date_breaks = "1 month", labels = label_first_mth_and_januarys)

ppb_mean_med_plot = ggplot(melt.d) +
  geom_line(aes(x = date, y = anchor_val, group = serial), alpha = 0.01) +
  geom_line(data = melt.d[, .(med = median(anchor_val), mean = mean(anchor_val)), by = .(date)],
            aes(x = date, y = med, color = "Median"), size = 2) +
  geom_line(data = melt.d[, .(mean = mean(anchor_val)), by = .(date)],
            aes(x = date, y = mean, color = "Mean"), size = 2, lty = "11") +
  theme_light() +
  theme(legend.position = 'bottom', text = element_text(size = 20)) +
  ylim(c(0,1)) +
  labs(x = "Date", y = "PPB value", color = "Central tendency metric") +
  scale_x_date(date_breaks = "1 month", labels = label_first_mth_and_januarys) +
  scale_color_manual(breaks = c("Mean", "Median"), values = c("lightblue", "darkred"))

ggsave(filename = file.path(output_path, "1k_ppb_fits.png"), plot = ppb_plot, device = 'png', units = 'in', height = 10, width = 15, dpi = 300)
ggsave(filename = file.path(output_path, "1k_ppb_fits_mean_med.png"), plot = ppb_mean_med_plot, device = 'png', units = 'in', height = 10, width = 15, dpi = 300)

rm(list=ls())

library(RSQLite)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(lubridate)

#### Initialization and Database ----
totalpop_10k <- 375.5 # Adjust accordingly
#totalpop_10k <- 2061.6 # Adjust accordingly
max_week <- 82 # Adjust accordingly, start counting from 1
con <- dbConnect(SQLite(), "./covid_vac_v7.0.sqlite")
pal <- c("#000000FF", "#440154FF", "#2C5FFFFF", 
         "#238A8DFF", "#49C16DFF", "#C8CB15FF")

dbListTables(con)
met <- dbGetQuery(con, "SELECT * FROM met")
#pars <- dbGetQuery(con, "SELECT * FROM par")

## Lookup table for week -> date
start_date <- ymd("2020-02-05")
dates <- seq(start_date, start_date + max_week*7, by = "week")
weeks_df <- data.frame(week = -1:(max_week-1), date = dates)
months <- dates
day(months) <- 1
months <- unique(months)

#### Functions ----
calc_eff <- function (incd, incd_control, start_wk) {
  # start_wk = when should the cumulation begin, range: 1 to Nweek
  Nweek <- length(incd)
  eff2 <- rep(NA, Nweek)
  if (start_wk > 1) eff2[1:(start_wk-1)] <- 0
  eff2[start_wk:Nweek] <- 
    (1 - cumsum(incd[start_wk:Nweek]) / cumsum(incd_control[start_wk:Nweek]))
  return(eff2)
}

plot_curves <- function (epc, param, legend = F, hline = NULL, ...) {
  
  flag <- ifelse(startsWith(param, "Eff") | param == "case_avert" |
                   param == "death_avert", 
                 "n", "l")

  cond <- epc$vac == 0
  plot(epc[cond, c("date", param)], 
       type = "n",
       col = pal[1],
       lwd = 2,
       axes = F, ...)
  axis(side = 1, 
       at = months, 
       labels = format(months, "%b-%y"),
       cex = 0.75)
  axis(side = 2)
  
  cond <- epc$vac == 1 & epc$vac_rate == 0
  lines(epc[cond, c("date", param)],
        col = pal[3],
        lwd = 2)
  
  cond <- epc$vac == 1 & epc$vac_rate == 1
  lines(epc[cond, c("date", param)],
        col = pal[5],
        lwd = 2)
  
  if (flag == "l") {
    cond <- epc$vac == 0
    lines(epc[cond, c("date", param)],
          col = pal[1],
          lwd = 2)
  }
  
  if (legend) {
    legend("right",
           legend = c("No Vaccine", "Slower rollout", 
                      "Faster rollout"),
           col = pal[c(1,3,5)],
           lty = 1,
           lwd = 2,
           bty = "n",
           xpd = NA)
  }
  
  if (!is.null(hline)) abline(h = hline, lty = 3, lwd = 2, xpd = F)
}

#### Calculate Effectiveness, Averted and Stuff
## Read from sqlite
sql <- paste0("SELECT vac, vac_rate, realization, mutation, serial FROM par")
par <- dbGetQuery(con, sql)

## Joining table
dat <- par %>%
  left_join(met, by = "serial")

#### Pivot longer ----
## Using vac == 0 as control for all other scenarios
epi_long <- dat %>%
  pivot_longer(-c(vac, vac_rate, mutation, realization, serial),
               names_to = c(".value", "week"),
               names_pattern = "([a-zA-Z]+)([0-9]+)")

epi_long <- epi_long %>%
  mutate(incd_c = c / totalpop_10k,
         incd_d = d / totalpop_10k,
         week = as.numeric(week))

#### Calculate effectiveness and case averted ----
### Eff1 = Cumulative from beginning
### Eff2 = Cumulative since vaccination
### We can actually do realization matching here

## Add incidence of control
epi_long1 <- epi_long %>%
  group_by(mutation, realization, week) %>%
  mutate(incd_c_control = incd_c[vac == 0],
         incd_d_control = incd_d[vac == 0]) %>%
  arrange(realization, week)

## For each trt combination and realization,
## calculate effectiveness & case averted
start_wk <- 47
epi_long1 <- epi_long1 %>%
  group_by(vac, vac_rate, mutation, realization) %>%
  mutate(Eff1_c = calc_eff(incd_c, incd_c_control, 1),
         Eff2_c = calc_eff(incd_c, incd_c_control, start_wk),
         Eff1_d = calc_eff(incd_d, incd_d_control, 1),
         Eff2_d = calc_eff(incd_d, incd_d_control, start_wk),
         case_avert = cumsum(incd_c_control) - cumsum(incd_c),
         death_avert = cumsum(incd_d_control) - cumsum(incd_d))

## Take median of all metrics along all realization
epi_long1 <- epi_long1 %>%
  group_by(vac, vac_rate, mutation, week) %>%
  summarise(incd_c = median(incd_c),
            incd_d = median(incd_d),
            Rt = median(Rt),
            Eff1_c = median(Eff1_c),
            Eff2_c = median(Eff2_c),
            Eff1_d = median(Eff1_d),
            Eff2_d = median(Eff2_d),
            case_avert = median(case_avert),
            death_avert = median(death_avert))

epi_long1 <- epi_long1 %>%
  left_join(weeks_df, by = "week")


for (m in 0:1) {
  epc <- epi_long1 %>%
    ungroup() %>%
    filter(mutation == m)
  
  outfile <- paste0("fig/fl-vac-mutation-", m, ".png")
  
  ## Plot start
  png(outfile, width = 20, height = 16, units = "cm",
      res = 200)
  
  par(mfrow = c(4, 2),
      mar = c(2, 4, 2, 1),
      oma = c(0, 0, 3, 0))
  
  plot_curves(epc, param = "Rt", legend = F, hline = 1, 
              xlab = "",
              ylab = "", 
              # xlim = c(0, max_week - 1),
              ylim = c(min(epi_long1$Rt), max(epi_long1$Rt)), 
              main = "Rt")
  plot.new()
  legend("center",
         legend = c("No Vaccine", "Slower rollout", 
                    "Faster rollout"),
         col = pal[c(1,3,5)],
         lty = 1,
         lwd = 2,
         bty = "n",
         xpd = NA,
         cex = 1.3)
  plot_curves(epc, param = "incd_c", legend = F,
              xlab = "",
              ylab = "", 
              # xlim = c(0, max_week - 1),
              ylim = c(0, max(epi_long1$incd_c)), 
              main = "Case incidence per 10k")
  plot_curves(epc, param = "incd_d", legend = F,
              xlab = "",
              ylab = "", 
              # xlim = c(0, max_week - 1),
              ylim = c(0, max(epi_long1$incd_d)), 
              main = "Death incidence per 10k")
  plot_curves(epc, param = "Eff2_c", legend = F, hline = 0, 
              xlab = "",
              ylab = "", 
              # xlim = c(0, max_week - 1),
              ylim = range(epi_long1$Eff2_c), 
              main = "Cumu. effectiveness (case)")
  plot_curves(epc, param = "Eff2_d", legend = F, hline = 0, 
              xlab = "",
              ylab = "", 
              # xlim = c(0, max_week - 1),
              ylim = range(epi_long1$Eff2_d), 
              main = "Cumu. effectiveness (death)")
  plot_curves(epc, param = "case_avert", legend = F, hline = 0, 
              xlab = "",
              ylab = "", 
              # xlim = c(0, max_week - 1),
              ylim = range(epi_long1$case_avert), 
              main = "Cumu. case averted per 10k")
  plot_curves(epc, param = "death_avert", legend = F, hline = 0, 
              xlab = "",
              ylab = "", 
              # xlim = c(0, max_week - 1),
              ylim = range(epi_long1$death_avert), 
              main = "Cumu. death averted per 10k")
  
  main <- ifelse(m == 0, "No mutant strain", "With mutant strain")
  
  # mtext("Week after first introduction", 1, 1, outer = T)
  mtext(main, 3, 0, outer = T)
  
  dev.off()
}

dbDisconnect(con)

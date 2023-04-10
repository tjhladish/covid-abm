library(data.table)
library(lubridate)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/state-vac/") }

#' ACS_2019_pop_data.csv available in active vac dir
#' cdc_covid-19_vax_data.csv available from the download_cdc_covid_vax_data.sh script
#' United_States_COVID-19_Cases_and_Deaths_by_State_over_Time.csv available from https://data.cdc.gov/Case-Surveillance/United-States-COVID-19-Cases-and-Deaths-by-State-o/9mfq-cb36
#' COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries.csv available from https://healthdata.gov/Hospital/COVID-19-Reported-Patient-Impact-and-Hospital-Capa/g62h-syeh/data 
.args <- if (interactive()) c(
  "ACS_2019_pop_data.csv",
  "cdc_covid-19_vax_data.csv",
  "United_States_COVID-19_Cases_and_Deaths_by_State_over_Time.csv",
  "COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries.csv"
) else commandArgs(trailingOnly = TRUE)

output_path = file.path('.', 'fig', 'state_vax_manuscript')
dir.create(output_path, recursive = TRUE)

#' locations that are included in the analysis
locs_of_interest <- c("FL", "VT", "MS")

#' read in population data, empirical vaccination data
pop.in <- fread(.args[1], header = TRUE)

cdc.in <- fread(.args[2])
cdc.in[, Date := as.Date(Date, "%m/%d/%Y") ]

dths.in <- fread(.args[3])[state %in% locs_of_interest]
dths.in[, date := mdy(submission_date)]

hosp.in <- fread(.args[4])[state %in% locs_of_interest]
hosp.in[, date := ymd(date)]

#' function to calculate the correct age bin populations
extract_bins <- function(SD) {
  SD[, 
     .(pop = c(
       `5_14` + `15_17` + `18_120`,
       `5_9` + (2/5)*`10_14`,
       `5_14` + `15_17` - `5_9`,
       `18_120` - `65_120`,
       `65_120`),
       bin_min = c(5,5,12,18,65),
       bin_max = c(120, 11, 17, 64, 120)
     )
  ]
}

#' generate long-form population data.table and calculate age bin population proportions
pop.dt <- pop.in[location %in% locs_of_interest, extract_bins(.SD), by=location]
pop.dt[, prop := pop/pop[(bin_min == 5) & (bin_max == 120)], by=location]

adj_dths <- dths.in[, .(date, state, tot_death)][
  pop.dt[bin_min == 5 & bin_max == 120, .(location, pop)], on = .(state==location)][
    , tot_dth_p10k := tot_death * (1e4/pop)]
adj_hosp <- hosp.in[is.finite(previous_day_admission_adult_covid_confirmed) | is.finite(previous_day_admission_pediatric_covid_confirmed), .(date, state, hosp_inc = previous_day_admission_adult_covid_confirmed + previous_day_admission_pediatric_covid_confirmed)][
  pop.dt[bin_min == 5 & bin_max == 120, .(location, pop)], on = .(state==location)][
    , hosp_inc_p10k := hosp_inc * (1e4/pop)][
      order(date)][
      , tot_hosp_inc_p10k := cumsum(hosp_inc_p10k), by = .(state)]

#' function to roll negative values over to the next day to address errors in cumulative data
rollneg <- function(v) {
  for (i in 1:(length(v) - 1)) {
    if (v[i] < 0) {
      n_to_decrement = v[i]
      v[i] = 0
      v[i+1] = v[i+1] + n_to_decrement
    }
  }
  return(as.numeric(v))
}

#' master function to process vaccination administration data
MASTERFUN <- function(SD) {
  #' subset CDC data to isolate necessary columns and fill NA values with 0
  sub_d = SD[Recip_Administered > 0 & Administered_Dose1_Recip_65Plus != 0,
             .(date = Date,
               Administered_Dose1_Recip_5Plus, Administered_Dose1_Recip_12Plus, Administered_Dose1_Recip_18Plus, Administered_Dose1_Recip_65Plus,
               Series_Complete_5Plus, Series_Complete_12Plus, Series_Complete_18Plus, Series_Complete_65Plus,
               Additional_Doses_12Plus, Additional_Doses_18Plus, Additional_Doses_65Plus)]
  sub_d = setnafill(x = sub_d, fill = 0, cols = colnames(sub_d)[2:ncol(sub_d)])
  
  #' because CDC data is reported using nested age bins, an adjustment is necessary for proper processing
  #' larger age bins may have 0s reported, but the following lines ensure that larger age bins have values that are at least equal to nested values
  sub_d[Administered_Dose1_Recip_18Plus == 0, Administered_Dose1_Recip_18Plus := Administered_Dose1_Recip_65Plus]
  sub_d[Administered_Dose1_Recip_12Plus == 0, Administered_Dose1_Recip_12Plus := Administered_Dose1_Recip_18Plus]
  sub_d[Administered_Dose1_Recip_5Plus == 0, Administered_Dose1_Recip_5Plus := Administered_Dose1_Recip_12Plus]
  
  sub_d[Series_Complete_18Plus == 0, Series_Complete_18Plus := Series_Complete_65Plus]
  sub_d[Series_Complete_12Plus == 0, Series_Complete_12Plus := Series_Complete_18Plus]
  sub_d[Series_Complete_5Plus == 0, Series_Complete_5Plus := Series_Complete_12Plus]
  
  sub_d[Additional_Doses_18Plus == 0, Additional_Doses_18Plus := Additional_Doses_65Plus]
  sub_d[Additional_Doses_12Plus == 0, Additional_Doses_12Plus := Additional_Doses_18Plus]
  
  #' this block address multiple processing steps
  #' nested to bounded age bins
  #' wide to long transformation
  #' cumulative to daily diffs (with the first value appended to ensure that total sum is maintained)
  #' apply rollneg function to ensure data is properly cumulative 
  diff.dt <- melt(sub_d[, .(
    date,
    cumul_65_120_dose_1 = Administered_Dose1_Recip_65Plus,
    cumul_18_64_dose_1 = Administered_Dose1_Recip_18Plus - Administered_Dose1_Recip_65Plus,
    cumul_12_17_dose_1 = Administered_Dose1_Recip_12Plus - Administered_Dose1_Recip_18Plus,
    cumul_5_11_dose_1 = Administered_Dose1_Recip_5Plus - Administered_Dose1_Recip_12Plus,
    
    cumul_65_120_dose_2 = Series_Complete_65Plus,
    cumul_18_64_dose_2 = Series_Complete_18Plus - Series_Complete_65Plus,
    cumul_12_17_dose_2 = Series_Complete_12Plus - Series_Complete_18Plus,
    cumul_5_11_dose_2 = Series_Complete_5Plus - Series_Complete_12Plus,
    
    cumul_65_120_dose_3 = Additional_Doses_65Plus,
    cumul_18_64_dose_3 = Additional_Doses_18Plus - Additional_Doses_65Plus,
    cumul_12_17_dose_3 = Additional_Doses_12Plus - Additional_Doses_18Plus
  )], id.vars = "date")[
    order(date), 
    .(date, value = c(value[1], diff(value))), by=.(variable = gsub("cumul_","daily_", variable))
  ][,
    .(date, value = rollneg(value)), by=variable
  ]
  return(diff.dt)
}

#' apply MASTERFUN to each location and clean up any remaining negative values at the end of the resulting TS
cdc.extract <- cdc.in[Location %in% locs_of_interest, MASTERFUN(.SD), by=.(location = Location)]
cdc.extract[date == max(date), value := ifelse(value < 0, 0, value)]

#' extract necessary columns from the variable column
cdc.extract[, c("bin_min", "bin_max", "dose") := tstrsplit(variable, "_", keep=c(2,3,5))]

#' join pop.dt to cdc.extract to calculate necessary proportions by empirical location
#' bin_prop: % of that age bin vaxd that day
#' tot_prop: % of total pop vaxd that day
dt <- cdc.extract[, .(location, date, bin_min=as.integer(bin_min), bin_max=as.integer(bin_max), dose=as.integer(dose), value)][
  order(location, date)][
    pop.dt, on=.(location, bin_min, bin_max)]
dt[, bin_prop := value / pop]
dt[, tot_prop := value / pop[bin_min == 5 & bin_max == 120], by=location]

dt <- dt[!is.na(value)]
dt[, c("pop", "prop") := NULL]

label_first_mth_and_januarys <- function(x) {
  ifelse(is.na(lag(x)) | year(lag(x)) != year(x),
         paste0(month(x, label = T), "\n", year(x)),
         paste0(month(x, label = T))
  )
}

shr <- list(
  theme_light(),
  xlab(""),
  theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), text = element_text(size = 15)),
  scale_color_manual(name = "", breaks = c("MS", "FL", "VT"), values = c("#F8766D", "#00BA38", "#619CFF")),
  scale_x_date(date_breaks = "1 month", labels = label_first_mth_and_januarys),
  ylim(c(0, NA))
)

cov = ggplot() +
  geom_line(data = dt[order(location, date, dose), .(tot_pop_vaxd = sum(tot_prop)), by = .(location, date, dose)][
    ,cumsum := cumsum(tot_pop_vaxd), by=.(location, dose)][
      ,.(date, cumsum), by = .(location, dose)][dose == 1],
    aes(x = date, y = cumsum, color = location), size = 1) +
  shr +
  labs(title = "Dose 1 Coverage", y = "") +
  theme(legend.position = c(0.99, 0.01), legend.justification = c("right", "bottom"), legend.direction = "vertical")
  

dths = ggplot(adj_dths[date >= ymd("2021-03-05") & date <= ymd("2022-08-03")]) +
  geom_line(aes(x = date, y = tot_dth_p10k, color = state), size = 1) +
  shr +
  labs(title = "Cumulative Deaths per 10K", y = "") +
  theme(legend.position = "none")

hosp = ggplot(adj_hosp[date >= ymd("2021-03-05") & date <= ymd("2022-08-03")]) +
  geom_line(aes(x = date, y = tot_hosp_inc_p10k, color = state), size = 1) +
  shr +
  labs(title = "Cumulative Hospitalizations per 10K", y = "") +
  theme(legend.position = "none")
  
final = cov / hosp / dths

ggsave(filename = file.path(output_path, "state_dth_vax_comp.png"), plot = final, device = png,
       units = "px", width = 3600, height = 3600, dpi = 360)

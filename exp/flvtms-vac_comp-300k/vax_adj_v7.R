library(data.table)
library(lubridate)
library(ggplot2)
library(viridis)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/flvtms-vac_comp-300k/") }

.args <- if (interactive()) c(
  "ACS_2019_pop_data.csv",
  "cdc_covid-19_vax_data.csv",
  "./dose_data"
) else commandArgs(trailingOnly = TRUE)

output_path = file.path('.', 'fig', 'state_vax_input')
dir.create(output_path, recursive = TRUE)

#' read in population data, empirical vaccination data
pop.in <- fread(.args[1], header = TRUE)
cdc.in <- fread(.args[2])
cdc.in[, Date := as.Date(Date, "%m/%d/%Y") ]

#' locations that are included in the analysis
locs_of_interest <- c("FL", "VT", "MS")

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

#' join pop.dt to dt to calculate similar proportions by simulation experiment
#' exp_value: multiply the state-specific age-bin vax proportion by FL's age-bin pop to calculate the adjuisted number of people vaxd in that bin
#' recalculate bin_prop and tot_prop using adjusted exp_value
exp.dt <- dt[pop.dt[location=="FL"], on=.(bin_min, bin_max)]
exp.dt[, `:=` (
  exp = paste0("FL_like_", location),
  exp_value = bin_prop * pop
)]
exp.dt[, exp_bin_prop := exp_value / pop]
exp.dt[, exp_tot_prop := exp_value / pop[bin_min == 5 & bin_max == 120]]

exp.dt <- exp.dt[!is.na(exp_value), .(exp, ref_loc = location, date, bin_min, bin_max, dose, value=exp_value, bin_prop=exp_bin_prop, tot_prop=exp_tot_prop)]
exp.dt[, value_p10k := bin_prop * 1e4]
exp.dt[, tot_doses_avail_p10k := sum(tot_prop * 1e4), by = .(exp, date)]

#' save subset of the transformed data to send to the ABM in terms of per 10k people
#' date, bin_min, bin_max, num_dose1_administered_p10k
q.dt <- exp.dt[dose == 1, .(date, ref_loc, bin_min, bin_max, bin_dose_1_p10k = (bin_prop * 1e4)), by = exp]

#' go from day 0 up to and including first data day in new dataset
#' divide by max in old set
#' multiply by first day in new set and remove first day
#' pre-pend new transformed series to new dataset to complete TS
dose.in = data.table()
for (loc in locs_of_interest) {
  tmp <- fread(paste0(.args[3], "/trends_in_number_of_covid19_vaccinations_in_", tolower(loc), ".csv"))
  dose_table <- tmp[, .(location = loc, date = Date, first = `Daily Count People Receiving Dose 1`,
                       second = `Daily Count of People Fully Vaccinated`, third = `Daily Count People Receiving a First Booster Dose`)]
  dose.in <- rbindlist(list(dose.in, dose_table))
}
dose.in[pop.dt[bin_min == 5 & bin_max == 120, .(location, pop)], on = .(location), tot_pop := pop][, tot_doses := sum(first+second+third), by = .(location, date)]
dose.in[, `:=` (
  `1_p10k` = first * (1e4/tot_pop),
  `2_p10k` = second * (1e4/tot_pop),
  `3_p10k` = third * (1e4/tot_pop)
)]

dose.in <- melt(dose.in[, .(location, date, `1_p10k`, `2_p10k`, `3_p10k`)], id = c("location", "date"), variable.name = "dose", value.name = "emp_value_p10k")
dose.in[, dose := tstrsplit(dose, "_", keep=1)][, dose := as.integer(dose)]

dose.dt <- exp.dt[, .(date, location = ref_loc, bin_min, bin_max, dose, value_p10k)][dose.in, on = .(location, date, dose)]

prepend_vax_doses <- function(SD, orig_ts) {
  tmp_ts <- cumsum(orig_ts)
  if (max(tmp_ts) == 0) {
    tmp_ts <- rep(0, length(tmp_ts))
  } else {
    tmp_ts <- (tmp_ts / max(tmp_ts)) * SD[date == min(date[!is.na(value_p10k)]), value_p10k]
  }
  tmp_ts <- c(tmp_ts[1], diff(tmp_ts), SD[date > min(date[!is.na(value_p10k)]), value_p10k])
  return(tmp_ts)
}

prepended_doses.dt = data.table()
for (loc in locs_of_interest) {
  for (d in 1:3) {
    orig_prepend_ts = dose.dt[location == loc & dose == d & is.na(value_p10k), emp_value_p10k]
    orig_prepend_ts = c(orig_prepend_ts, dose.dt[location == loc & dose == d & date == min(date[!is.na(value_p10k)]), unique(emp_value_p10k)])
    for (bin in unique(dose.dt[!is.na(bin_min), bin_min])) {
      if (bin == 5 & d == 3) { next }
      tmp = dose.dt[location == loc & dose == d & bin_min == bin][
        data.table(date = seq.Date(as.Date("2020-12-14"), as.Date("2022-03-07"), by = "day")), on = "date"]
      tmp[, `:=` (
        location = unique(location[!is.na(location)]),
        bin_min = unique(bin_min[!is.na(bin_min)]),
        bin_max = unique(bin_max[!is.na(bin_max)]),
        dose = unique(dose[!is.na(dose)])
      )]
      tmp[, adj_value_p10k := prepend_vax_doses(.SD, orig_prepend_ts)]
      prepended_doses.dt = rbindlist(list(prepended_doses.dt, tmp[, .(date, location, bin_min, bin_max, dose, n_doses_p10k=adj_value_p10k)]))
    }
  }
}

#' output queue file and dosing file to current directory with no urgent doses allocated
out.dt <- prepended_doses.dt[,.(date, ref_location=location, bin_min, bin_max, dose, n_doses_p10k)]
out.dt[, is_urg := 0]
tmp <- prepended_doses.dt[,.(date, ref_location=location, bin_min, bin_max, dose, n_doses_p10k)]
tmp[, `:=` (
  is_urg = 1,
  n_doses_p10k = 0
)]
out.dt <- rbindlist(list(out.dt, tmp))
fwrite(out.dt[,.(date, ref_location, bin_min, bin_max, dose, is_urg, n_doses_p10k)], file = "state_based_counterfactual_doses.txt", sep = ' ')

active_dosing = prepended_doses.dt[,.(date, ref_location=location, bin_min, bin_max, dose, n_doses_p10k)][ref_location == "FL"]
active_dosing[, is_urg := 0]

#' the chunk below will allocate the same amount of first doses used empirically after May 1, 2021 for the active campaign
# tmp <- prepended_doses.dt[,.(date, ref_location=location, bin_min, bin_max, dose, n_doses_p10k)][ref_location == "FL"]
# tmp[, is_urg := 1]
# tmp[dose != 1 | date < "2021-05-01", n_doses_p10k := 0]

#' the chunk below will allocate a constant amount of first doses for the active campaign after May 1, 2021
#const_allocation = 5e4 * (1e4/pop.in[location == "FL", `0_17`+`18_120`])
const_allocation = 1e5 * (1e4/pop.in[location == "FL", `0_17`+`18_120`])
tmp <- prepended_doses.dt[,.(date, ref_location=location, bin_min, bin_max, dose, n_doses_p10k)][ref_location == "FL"]
tmp[, `:=` (is_urg = 1, n_doses_p10k = 0)]
tmp[dose == 1 & bin_min == 5, n_doses_p10k := const_allocation]

#' the chunk below will allocate XX% of VT-scenario doses to the passive and 1-XX% to the active to ensure that the passive reaches the emp coverage endpoint in FL
# tmp <- prepended_doses.dt[,.(date, ref_location=location, bin_min, bin_max, dose, n_doses_p10k)][ref_location == "VT" | ref_location == "FL"]
# tmp_cov <- tmp[, .(tot_doses = sum(n_doses_p10k)), by =.(ref_location, dose, bin_min)]
# adj <- tmp_cov[, .(passive_prop = tot_doses[ref_location == "FL"]/tot_doses[ref_location == "VT"]), by = .(dose, bin_min)][, active_prop := 1 - passive_prop]
# tmp <- tmp[ref_location == "VT"][adj, on = .(dose, bin_min)]
# tmp[, `:=` (passive_doses = n_doses_p10k * passive_prop, active_doses = n_doses_p10k * active_prop, is_urg = 1)]
# active_dosing <- tmp[, .(date, ref_location, bin_min, bin_max, dose, n_doses_p10k = passive_doses, is_urg)][, `:=` (ref_location = "FL", is_urg = 0)]
# tmp <- tmp[, .(date, ref_location, bin_min, bin_max, dose, n_doses_p10k = active_doses, is_urg)][, ref_location := "FL"]

active_dosing <- rbindlist(list(active_dosing, tmp))
fwrite(active_dosing[,.(date, ref_location, bin_min, bin_max, dose, is_urg, n_doses_p10k)], file = "active_vax_counterfactual_doses.txt", sep = ' ')

# ggplot() +
#   geom_line(data = prepended_doses.dt[, .(date, cumsum=cumsum(n_doses_p10k)), by=.(location, bin_min, dose)], aes(x = date, y = cumsum, color=location), linetype='dashed') +
#   geom_line(data = exp.dt[, .(date, cumsum=cumsum(value_p10k)), by=.(location=ref_loc, bin_min, dose)], aes(x = date, y = cumsum, color=location)) +
#   geom_line(data = active_dosing[, .(doses=sum(n_doses_p10k)), by = .(date, dose, bin_min, is_urg)][, cumsum := cumsum(doses), by = .(dose, bin_min, is_urg)][is_urg == 0],
#             aes(x = date, y = cumsum), size = 0.5, lty = "dashed", color = "red") +
#   geom_line(data = active_dosing[, .(doses=sum(n_doses_p10k)), by = .(date, dose, bin_min)][, cumsum := cumsum(doses), by = .(dose, bin_min)],
#             aes(x = date, y = cumsum), size = 0.5, color = "black", lty = "dashed") +
#   geom_hline(yintercept = 1e4, linetype = 'dotted') +
#   facet_grid(rows=vars(dose), cols=vars(bin_min)) +
#   xlim(c(as.Date("2020-12-14"), as.Date("2022-03-01")))

#' fancy stuff from CABP
#' gg.scale.wrapper <- function(
#' 	scale_fun,
#' 	...
#' ) {
#' 	stopifnot(!missing(scale_fun))
#' 	defs <- list(...)
#' 	if (!length(defs)) warning(
#' 		"provided no default arguments; consider using scale_fun directly."
#' 	)
#' 
#' 	return(function(...) {
#' 		#' this different ... is how you get a function back that let's you
#' 		#' override defaults, set other arguments to scale_... functions
#' 		.ellipsis <- list(...)
#' 		.args <- defs
#' 		.args[names(.ellipsis)] <- .ellipsis
#' 		do.call(scale_fun, .args)
#' 	})
#' }
#' 
#' scale_color_location <- gg.scale.wrapper(
#' 	scale_color_viridis_d,
#' 	name = "Location", end = 0.9
#' )

shr <- list(
	scale_x_date(date_breaks = "months", date_labels = "%b"),
	xlab("Date (2021-2022)"),
	ylab("Vaccination coverage"),
	theme_light(),
	theme(legend.position = 'bottom', panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())
)

adj_binned_cov_fig <- ggplot() +
  geom_line(data = exp.dt[order(exp, date, dose), .(tot_bin_vaxd = sum(bin_prop)), by = .(exp, date, dose, bin_min)][
                          ,cumsum := cumsum(tot_bin_vaxd), by=.(exp, dose, bin_min)][
                            ,.(date, cumsum), by = .(exp, dose, bin_min)],
            aes(x = date, y = cumsum, color = as.factor(exp))) +
  geom_line(data = dt[order(location, date, dose), .(tot_bin_vaxd = sum(bin_prop)), by = .(location, date, dose, bin_min)][
                        ,cumsum := cumsum(tot_bin_vaxd), by=.(location, dose, bin_min)][
                          ,.(date, cumsum), by = .(location, dose, bin_min)],
            aes(x = date, y = cumsum, color = location), linetype = 'dashed') +
  geom_hline(yintercept = 1, linetype = 'dotted') +
  scale_color_manual(name = "Experiment",
                     breaks = c("FL_like_FL", "FL_like_VT", "FL_like_MS", "FL", "MS", "VT"),
                     labels = c("Baseline", "FL like VT", "FL like MS", "FL (emp.)", "MS (emp.)", "VT (emp.)"),
                     values = c(viridis(3, end = 0.9), 'black', 'black', 'black')) +
  facet_grid(rows = vars(dose), cols = vars(as.integer(bin_min))) +
  ggtitle("Total coverage per age bin per dose (empirical vs. adjusted)") +
  shr

tot_cov_comparison <- ggplot() +
  geom_line(data = dt[order(location, date, dose), .(tot_pop_vaxd = sum(tot_prop)), by = .(location, date, dose)][
                      ,cumsum := cumsum(tot_pop_vaxd), by=.(location, dose)][
                        ,.(date, cumsum), by = .(location, dose)],
            aes(x = date, y = cumsum)) +
  geom_line(data = dt[order(location, date, dose), .(tot_pop_vaxd = sum(tot_prop)), by = .(location, date, dose)][
                      ,cumsum := cumsum(tot_pop_vaxd), by=.(location, dose)][
                        location == "FL", .(date, cumsum), by = .(dose)],
            aes(x = date, y = cumsum), color=viridis(3, end = 0.8)[3]) +
  geom_line(data = exp.dt[order(exp, date, dose), .(tot_pop_vaxd = sum(tot_prop)), by = .(exp, date, dose)][
                          ,cumsum := cumsum(tot_pop_vaxd), by=.(exp, dose)][
                            ,.(date, location = sub("FL_like_(.*)", "\\1", exp), cumsum), by = .(exp, dose)],
            aes(x = date, y = cumsum), color=viridis(3, end = 0.8)[3], linetype='dashed') +
  geom_hline(yintercept = 1, linetype = 'dotted') +
  facet_grid(rows = vars(dose), cols = vars(location)) +
  ggtitle("Total coverage per state per dose (empirical vs. adjusted)") +
  shr

bin.labs <- c("5-11", "12-17", "18-64", "65+")
names(bin.labs) <- c(5, 12, 18, 65)

dose_delivery <- ggplot() +
  geom_line(data = prepended_doses.dt[, .(date, cumsum=cumsum(n_doses_p10k)), by=.(location, bin_min, dose)], aes(x = date, y = cumsum, color=location), linetype='dashed') +
  geom_line(data = exp.dt[, .(date, cumsum=cumsum(value_p10k)), by=.(location=ref_loc, bin_min, dose)], aes(x = date, y = cumsum, color=location)) +
  geom_hline(yintercept = 1e4, linetype = 'dotted') +
  facet_grid(rows=vars(dose), cols=vars(bin_min), labeller = labeller(bin_min = bin.labs)) +
  shr +
  scale_color_discrete(name = "Reference location") +
  xlab("Date (2020-2022)") + 
  ylab("Total doses delivered per 10K") +
  ggtitle("Cumulative administered doses (with prepending adjustemt)")

# ggplot() +
#   geom_line(data = prepended_doses.dt[, .(sum=sum(n_doses_p10k)), by=.(date, location, dose)][,.(date, cumsum=cumsum(sum)), by=.(location, dose)],
#             aes(x = date, y = cumsum)) +
#   geom_line(data = dose.in[, .(date, cumsum=cumsum(emp_value_p10k)), by=.(location, dose)], aes(x = date, y = cumsum), color='red') +
#   geom_line(data = exp.dt[order(exp, date, dose), .(tot_pop_vaxd = sum(value_p10k)), by = .(exp, date, dose)][
#     ,cumsum := cumsum(tot_pop_vaxd), by=.(exp, dose)][
#       ,.(date, location = sub("FL_like_(.*)", "\\1", exp), cumsum), by = .(exp, dose)],
#     aes(x = date, y = cumsum), color=viridis(3, end = 0.8)[3], linetype='dashed') +
#   facet_grid(rows = vars(dose), cols = vars(location))
# 
# ggplot() +
#   geom_line(data = exp.dt[order(exp, date, dose), .(tot_bin_vaxd = sum(bin_prop)), by = .(exp, date, dose, bin_min)][
#     ,cumsum := cumsum(tot_bin_vaxd), by=.(exp, dose, bin_min)][
#       ,.(date, cumsum), by = .(exp, dose, bin_min)][exp=="FL_like_FL"],
#     aes(x = date, y = cumsum, color = as.factor(exp))) +
#   geom_line(data = dt[order(location, date, dose), .(tot_bin_vaxd = sum(bin_prop)), by = .(location, date, dose, bin_min)][
#     ,cumsum := cumsum(tot_bin_vaxd), by=.(location, dose, bin_min)][
#       ,.(date, cumsum), by = .(location, dose, bin_min)][location=="FL"],
#     aes(x = date, y = cumsum, color = location), linetype = 'dashed') +
#   geom_line(data = dose.in[, .(date, cumsum=cumsum(emp_value_p10k/1e4)), by=.(location, dose)][location=="FL"], aes(x = date, y = cumsum), color='red') +
#   geom_line(data = prepended_doses.dt[, .(date, cumsum=cumsum(n_doses_p10k/1e4)), by=.(location, bin_min, dose)][location=="FL"],
#             aes(x = date, y = cumsum, color=location), linetype='dashed') +
#   scale_color_manual(name = "Experiment",
#                      breaks = c("FL_like_FL", "FL_like_VT", "FL_like_MS", "FL", "MS", "VT"),
#                      labels = c("Baseline", "FL like VT", "FL like MS", "FL (emp.)", "MS (emp.)", "VT (emp.)"),
#                      values = c(viridis(3, end = 0.9), 'black', 'black', 'black')) +
#   facet_grid(rows = vars(dose), cols = vars(as.integer(bin_min))) +
#   ggtitle("Total coverage per age bin per dose (empirical vs. adjusted)") +
#   shr

ggsave(filename = file.path(output_path, "tot_cov_comparison_v2.png"),        plot = tot_cov_comparison,
    device = 'png', units = 'in', height = 6, width = 15, dpi = 300)
ggsave(filename = file.path(output_path, "adj_binned_cov_comparison_v2.png"), plot = adj_binned_cov_fig,
    device = 'png', units = 'in', height = 6, width = 15, dpi = 300)
ggsave(filename = file.path(output_path, "dosing_adjustment_v2.png"),         plot = dose_delivery,
    device = 'png', units = 'in', height = 6, width = 15, dpi = 300)

startup_check <- function() {
  req_pkg <- c("dplyr", "data.table", "argparser", "lubridate",
               "ggplot2", "tidyr")
  reqs <- sapply(req_pkg, require, quietly=T, character.only=T)
  
  if(!all(reqs)) {
    ms_pkg <- paste(req_pkg[!reqs], collapse = ",")
    stop("Cannot find the required packages: ", ms_pkg)
  }
}

parser <- function(argv = NULL) {
  p <- arg_parser("Plot CFR from input csv vs empirical data")
  p <- add_argument(p, "input_csv", help = "path to input csv")
  p <- add_argument(p, "--output", help = "path to output (png)",
                    default = "output.png")
  p <- add_argument(p, "--ini_weeks_to_remove", 
                    help = "number of weeks to remove from initial time series",
                    default = 0)
  p <- add_argument(p, "--emp_csv", 
                    help = "path to emipirical csv",
                    default = "fl_weekly_ll_derived_chd.csv")
  
  if (is.null(argv)) {
    argv <- parse_args(p)
  } else {
    argv <- parse_args(p, argv = argv)
  }
  
  return(argv)
}

load_input <- function(input_csv) {
  df <- fread(input_csv) %>%
    mutate(week = epiweek(date),
           week = ifelse(week != 53 & year(date) != 2020, 53+week, week))
  weekly <- df %>%
    group_by(week) %>%
    summarise(death_in = sum(deaths), case_in = sum(cases))
  
  return(weekly)
}

comb_input_emp <- function(input_csv, emp_csv) {
  input_weekly <- load_input(input_csv)
  emp_weekly <- fread(emp_csv)
  
  weeklies <- emp_weekly %>%
    left_join(input_weekly, by = "week")
  return(weeklies)
}

calc_cfr_long <- function(weeklies, ini_week_to_rm = 0) {
  weeklies <- tail(weeklies, nrow(weeklies) - ini_week_to_rm)
  weeklies_long <- weeklies %>%
    mutate(emp = death_n/n, inp = death_in/case_in) %>%
    select(week, week_end_date, emp, inp) %>%
    pivot_longer(c(emp, inp), names_to = "type", values_to = "cfr")
  return(weeklies_long)
}

#### Main
startup_check()
# Debug 
# args <- parser(c("cfr.csv", "-i", "2"))
args <- parser()

df <- comb_input_emp(args$input_csv, args$emp_csv)
df_long <- calc_cfr_long(df, args$ini_weeks_to_remove)
p1 <- ggplot(df_long) +
  geom_line(aes(x=week_end_date, y=cfr, colour=type)) +
  scale_colour_manual(name="", label=c("Empirical", "Modelled"),
                      values = c("black", "red")) +
  labs(x="", y="CFR", title="Empirical vs Modelled CFR") +
  theme_bw() +
  theme(legend.position = c(.99, .99),
        legend.justification = c(1, 1),
        legend.title = element_blank())

ggsave(args$output, p1, width = 15, height = 10, units = "cm")

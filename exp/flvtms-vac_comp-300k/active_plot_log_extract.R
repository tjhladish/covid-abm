library(data.table)
library(DBI)

if (interactive()) { setwd("~/Documents/work/covid-projs/active-vac-exp/")}

.args = if (interactive()) c(
  "active_vac_raw_output.tgz",
  "test.sqlite"
) else commandArgs(trailingOnly = TRUE)

if (length(.args) != 2) {
  stop("Rscript active_plot_log_extract.R [path to plotlog files] [path to output SQLite DB]")
}

input_log_dir = .args[1]
out_db_path = .args[2]

out.d = data.table()

print("Extracting data from plot_log files...")
pb = txtProgressBar(min = 0, max = length(list.files(input_log_dir)), initial = 0, style = 3)

cols_to_extract = c(
  "date",
  "inf",
  "symp_infs",
  "sevr_infs",
  "crit_infs",
  "all_deaths",
  "std_doses",
  "urg_doses"
)

for (i in 1:length(list.files(input_log_dir))) {
  in_plot_log = tolower(list.files(input_log_dir)[i])
  in_serial = as.integer(sub("plot_log(.*).csv", "\\1", in_plot_log))
  
  if (!file.exists(file.path(input_log_dir, in_plot_log))) { next }
  
  d = fread(file.path(input_log_dir, in_plot_log))
  
  sub.d = d[, ..cols_to_extract]
  sub.d[, serial := in_serial]
  sub.d = sub.d[, .(serial, date, inf, symp_infs, sevr_infs, crit_infs, all_deaths, std_doses, urg_doses)]
  
  out.d = rbindlist(list(out.d, sub.d))
  setTxtProgressBar(pb, i)
}
close(pb)

out.d[, date := as.character(date)]
# CREATE TABLE meta ( serial int not null, date text not null, inf real, symp real, sev real, crit real, deaths real, std_doses real, urg_doses real, primary key(serial, date));
setnames(out.d, c('symp_infs', 'sevr_infs', 'crit_infs', 'all_deaths'), c('symp', 'sev', 'crit', 'deaths'))

print("Writing data to new database...")
db = dbConnect(RSQLite::SQLite(), out_db_path)
dbWriteTable(db, "meta", out.d, append=T)
dbDisconnect(db)

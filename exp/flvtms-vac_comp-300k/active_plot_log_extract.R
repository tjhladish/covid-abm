library(data.table)
library(DBI)

if (interactive()) { setwd("~/Documents/work/covid-projs/active-vac-exp/")}

.args = if (interactive()) c(
  "active_vac_raw_output.tgz",
  "test.sqlite"
) else commandArgs(trailingOnly = TRUE)

if (length(.args) != 2) {
  stop("Rscript active_plot_log_extract.R [path to zipped archive] [path to output SQLite DB]")
}

out_db_path = .args[2]

print("Extracting plot_log files from zipped archive...")
tmp_dir = paste0("./tmp", paste0(sample(c(sample(LETTERS, 10, replace = T), sample(0:9, 10, replace = T))), collapse = ''))
dir.create(tmp_dir)
untar(.args[1], exdir = tmp_dir)

out.d = data.table()

print("Extracting data from plot_log files...")
pb = txtProgressBar(min = 0, max = length(list.files(tmp_dir)), initial = 0, style = 3)

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

for (i in 1:length(list.files(tmp_dir))) {
  in_plot_log = tolower(list.files(tmp_dir)[i])
  in_serial = as.integer(sub("plot_log(.*).csv", "\\1", in_plot_log))
  
  if (!file.exists(file.path(tmp_dir, in_plot_log))) { next }
  
  d = fread(file.path(tmp_dir, in_plot_log))
  
  sub.d = d[, ..cols_to_extract]
  sub.d[, serial := in_serial]
  sub.d = sub.d[, .(serial, date, inf, symp_infs, sevr_infs, crit_infs, all_deaths, std_doses, urg_doses)]
  
  out.d = rbindlist(list(out.d, sub.d))
  setTxtProgressBar(pb, i)
}
close(pb)

out.d[, date := as.character(date)]

print("Writing data to new database...")
db = dbConnect(RSQLite::SQLite(), out_db_path)
dbWriteTable(db, "meta", out.d)
dbDisconnect(db)
unlink(tmp_dir, recursive = T)

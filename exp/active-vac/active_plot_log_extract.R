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

# CREATE TABLE meta ( serial int not null, date text not null, inf real, symp real, sev real, crit real, deaths real, std_doses real, urg_doses real, primary key(serial, date));
db = dbConnect(RSQLite::SQLite(), out_db_path)
for (i in 1:length(list.files(input_log_dir))) {
    in_plot_log = tolower(list.files(input_log_dir)[i])
    in_serial = as.integer(sub("plot_log(.*).csv", "\\1", in_plot_log))

    if (!file.exists(file.path(input_log_dir, in_plot_log))) { next }

    sub.d = fread(file.path(input_log_dir, in_plot_log), select = cols_to_extract)
    sub.d = sub.d[, .(serial = in_serial, date, inf, symp_infs, sevr_infs, crit_infs, all_deaths, std_doses, urg_doses)]

    sub.d[, date := as.character(date)]
    setnames(sub.d, c('symp_infs', 'sevr_infs', 'crit_infs', 'all_deaths'), c('symp', 'sev', 'crit', 'deaths'))

    dbWriteTable(db, "meta", sub.d, append=T)

    setTxtProgressBar(pb, i)
}
close(pb)

dbDisconnect(db)

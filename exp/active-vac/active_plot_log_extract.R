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

plotlogs <- list.files(input_log_dir, full.names = TRUE)
serials <- as.integer(gsub("^.+plot_log(.*)\\.csv$", "\\1", plotlogs))

print("Extracting data from plot_log files...")
pb = txtProgressBar(min = 0, max = length(plotlogs), initial = 0, style = 3)

cols_to_extract = c(
  "date",
  "inf",
#  "symp_infs",
  "sevr_infs",
#  "crit_infs",
  "all_deaths",
  "std_doses",
  "urg_doses",
  "Rt"
)

# CREATE TABLE meta ( serial int not null, date text not null, inf real, symp real, sev real, crit real, deaths real, std_doses real, urg_doses real, primary key(serial, date));
db = dbConnect(RSQLite::SQLite(), out_db_path)
for (i in seq_along(plotlogs)) {
    in_plot_log = plotlogs[i]
    in_serial = serials[i]

    #' if (!file.exists(file.path(input_log_dir, in_plot_log))) { next }

    sub.d = fread(in_plot_log, select = cols_to_extract)
    #sub.d = sub.d[, .(serial = in_serial, date, inf, symp_infs, sevr_infs, crit_infs, all_deaths, std_doses, urg_doses)]
    sub.d = sub.d[, .(serial = in_serial, date, inf, sevr_infs, all_deaths, std_doses, urg_doses, Rt)]

    sub.d[, date := as.character(date)]
    #setnames(sub.d, c('symp_infs', 'sevr_infs', 'crit_infs', 'all_deaths'), c('symp', 'sev', 'crit', 'deaths'))
    setnames(sub.d, c('sevr_infs', 'all_deaths'), c('sev', 'deaths'))

    dbWriteTable(db, "meta", sub.d, append=T)

    setTxtProgressBar(pb, i)
}
close(pb)

dbDisconnect(db)

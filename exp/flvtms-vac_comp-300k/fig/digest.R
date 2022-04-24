
.pkgs <- c("data.table", "RSQLite", "jsonlite")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

.args <- if (interactive()) c(
  file.path("tmp", c(
    "abc_covid_counterfactuals.json", "covid_counterfactuals_v1.0.sqlite"
  )),
  file.path("fig", "digest.rds")
) else commandArgs(trailingOnly = TRUE)

js <- read_json(.args[1])

db <- dbConnect(SQLite(), .args[2])
stmt <- paste0("SELECT P.*, M.* FROM par AS P JOIN met AS M USING(serial);")
dt <- as.data.table(dbGetQuery(db, stmt))
dbDisconnect(db)

long.dt <- melt(
  dt[, .SD, .SDcols = unique(names(dt))],
  id.vars = c("serial", "seed", "vac", "realization", "mutation", "state")
)[,
  week := as.integer(gsub("[^0-9]+(\\d+)", "\\1", variable ))
][,
  variable := gsub("([^0-9]+)\\d+", "\\1", variable)
][variable %in% c("c", "s", "d")]

long.dt[, state := c("FL", "VT", "MS")[state+1] ]

saveRDS(long.dt, tail(.args, 1))

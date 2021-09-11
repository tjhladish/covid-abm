require(RSQLite)
require(data.table)

.args <- if (interactive()) c(
  "../../covid_vac_v9.0.sqlite",
  "rawresults.rds"
) else commandArgs(trailingOnly = TRUE)

con <- dbConnect(SQLite(), .args[1])
dt <- data.table(dbGetQuery(con, "SELECT M.*, P.* FROM met AS M JOIN par AS P USING(serial) JOIN job AS J USING(serial) WHERE J.status == 'D';"))
dbDisconnect(con)
# -1 drops duplicate serial column from grabbing both M.*, P.*
mlt <- melt(dt[, -1], id.vars = c("serial", "realization", "seed", "vac", "mutation"))
mlt[, week := as.integer(gsub("[^0-9]*([0-9]+)$","\\1",variable)) ]
mlt[, var := factor(gsub("([^0-9]*)[0-9]+$","\\1",variable)) ]

start_date <- as.Date("2020-02-05")

mlt[, date := start_date + week*7 ]

saveRDS(mlt, tail(.args, 1))

#' @examples 
#' require(ggplot2)
#' ggplot(mlt) + aes(week, value,
#'   color=c("WT","VAR")[mutation+1], linetype=c("NONE","VAC")[vac+1],
#'   group = interaction(serial,realization,vac,mutation)
#' ) +
#'   facet_grid(var ~ ., scales = "free") +
#'   geom_line(alpha = 0.1) +
#'   theme_minimal()

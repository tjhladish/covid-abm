
package <- c("data.table", "RSQLite")
pkg_cond <- sapply(package, require, character.only = T, quietly = T)
if (any(!pkg_cond)) stop("Following package(s) are required but not installed: ",
                         paste(package[!pkg_cond], collapse = ", "))
options(scipen = 10)

.args <- if (interactive()) c(
  "sim_pop/pop.sqlite",
  "sim_pop/locations.txt", "sim_pop/population.txt", "sim_pop/network.txt", "sim_pop/public-locations.txt"
) else commandArgs(trailingOnly = TRUE)

dbpth <- .args[1]

#### Connect
con <- DBI::dbConnect(RSQLite::SQLite(), dbpth)
extracurricular_mode <- "extracurr" %in% dbListTables(con)

#### Pulling query using SQL
## Locations
loc <- dbGetQuery(
  con, sprintf("
SELECT
  locid, ROUND(x, 5) AS x, ROUND(y, 5) AS y, type,
  COALESCE(nh.hf_locid, hf_locid2) - 1 AS hfid,
  ROUND(compliance, 4) AS compliance, wp.essential%s
FROM loc
LEFT JOIN nh USING (locid)
LEFT JOIN (
  SELECT locid, hf_locid AS hf_locid2, compliance FROM hh WHERE nh = 'n'
) USING (locid)
LEFT JOIN wp USING (locid);",
  ifelse(extracurricular_mode, ", transmission", "")
)) |> setDT()

loc[, locid := locid - 1]
loc[is.na(hfid), hfid := -1]
loc[type != "w", essential := "y"]
loc[is.na(compliance), compliance := -1]

if (extracurricular_mode) {
  loc[is.na(transmission), transmission := "N"]
}

message("Preview for loc")
print(head(loc))
fwrite(loc, tail(.args, 4)[1], sep = " ", scipen = 10)
rm(loc)

## Persons
pers <- dbGetQuery(con, "SELECT
pid, r.locid - 1 AS home_id,
sex, age, m.locid - 1 AS day_id
FROM pers
LEFT JOIN movement AS m
USING (pid)
LEFT JOIN reside AS r
USING (pid);") |> setDT()

pers[, pid := pid - 1]
pers[is.na(day_id), day_id := -1]

message("Preview for pers")
print(head(pers))
fwrite(pers, tail(.args, 4)[2], sep = " ", scipen = 10)
rm(pers)

## Household network (pretty much just take from the sql db)
hh_network <- dbGetQuery(
  con,
  "SELECT locid1 - 1 AS locid1, locid2 - 1 AS locid2 FROM hh_network;"
) |> setDT()

message("Preview for hh_network")
print(head(hh_network))
fwrite(hh_network, tail(.args, 4)[3], sep = " ", col.names = F, scipen = 10)
rm(hh_network)

## Extracurricular activities (pretty much just take from the sql db)
if (extracurricular_mode) {
  ec <- dbGetQuery(
    con,
    "SELECT dest_locid_1, dest_locid_2, dest_locid_3, dest_locid_4, dest_locid_5 FROM extracurr;"
  ) |> as.matrix() |> as.vector() |> sort() |> unique()
  ec <- ec - 1

  ec_df <- data.frame(ec = ec)
  message("Preview for public-locations")
  print(head(ec_df))
  fwrite(ec_df, tail(.args, 4)[4], sep = " ", col.names = F)
}

dbDisconnect(con)

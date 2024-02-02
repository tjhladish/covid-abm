rm(list=ls())

package <- c("dplyr", "data.table", "RSQLite", "stringr")
pkg_cond <- sapply(package, require, character.only = T, quietly = T)
if (any(!pkg_cond)) stop("Following package(s) are required but not installed: ",
                         paste(package[!pkg_cond], collapse = ", "))
options(scipen = 10)

#### Getting Rscript args
arg <- commandArgs(trailingOnly=TRUE)
arg <- trimws(arg)

if (length(arg) < 1) {
  stop("Argument (target tgz file name) is required.")
} else if (length(arg) > 1) {
  stop("Too many argument.")
} else if (str_detect(arg, "/")) {
  stop("Target tgz file should be in the same folder of this script, and you should execute this script
       in its directory, so argument with '/' character is inappropriate.")
} else if (!str_detect(arg, "sim_pop-.*-3.11.tgz")) {
  stop("Argument is not a valid sim_pop 3.11 tgz file.")
} else {
  tgz_file <- arg
}

#### Splicing args
cnt <- str_remove(tgz_file, "sim_pop-") %>% str_remove("-3.11.tgz")
folder <- str_remove(tgz_file, ".tgz")
sqlite_file <- str_replace(tgz_file, ".tgz", ".sqlite")
out_folder <- paste0("sim_pop-", cnt, "/")
dir.create(out_folder)

#### Unpack
message("Unpacking files...")
if (!file.exists(tgz_file)) stop("tgz file must be in the same folder as this script.")
untar(tgz_file, exdir=".")

#### Connect
con <- DBI::dbConnect(RSQLite::SQLite(), paste0(folder, "/", sqlite_file))
extracurricular_mode <- ifelse("extracurr" %in% dbListTables(con), T, F)

#### Pulling query using SQL
## Locations
sql <- "SELECT loc.locid, loc.x, loc.y, loc.type, nh.hf_locid FROM loc
LEFT JOIN nh
ON loc.locid = nh.locid"
loc <- dbGetQuery(con, sql)

hh_nonnh <- dbGetQuery(con, "SELECT locid, hf_locid AS hf_locid2, compliance FROM hh WHERE hh.nh = 'n'")
loc <- loc %>%
  left_join(hh_nonnh)
loc <- loc %>%
  mutate(hf_locid = ifelse(is.na(hf_locid), hf_locid2, hf_locid)) %>%
  select(-hf_locid2)

if (extracurricular_mode) {
  sql <- "SELECT locid, essential, transmission FROM wp"
  wp <- dbGetQuery(con, sql)
  addtl_col <- "transmission"
} else {
  sql <- "SELECT locid, essential FROM wp"
  wp <- dbGetQuery(con, sql)
  addtl_col <- "dummy"
}

loc <- loc %>%
  left_join(wp)
if ("hf_locid" %in% colnames(loc)) colnames(loc)[colnames(loc) == "hf_locid"] <- "hfid"
loc$essential[loc$type != "w"] <- "y"
loc$locid <- loc$locid - 1
loc$x <- round(loc$x, 5)
loc$y <- round(loc$y, 5)
loc <- loc %>%
  select(locid, x, y, type, essential, hfid, compliance, matches(addtl_col))

loc$hfid <- loc$hfid - 1
loc$hfid[is.na(loc$hfid)] <- -1
loc$compliance[is.na(loc$compliance)] <- -1
loc$compliance <- round(loc$compliance, 4)
if (extracurricular_mode) loc$transmission[is.na(loc$transmission)] <- "N"

message("Preview for loc")
print(head(loc))
fwrite(loc, paste0(out_folder, "locations-", cnt, ".txt"), sep = " ", scipen = 10)

## Persons
sql <- "SELECT p.pid, r.locid AS home_id, p.sex, p.age, m.locid AS day_id, p.undlycond
FROM pers AS p
LEFT JOIN movement AS m
ON p.pid = m.pid
LEFT JOIN reside AS r
ON p.pid = r.pid"
pers <- dbGetQuery(con, sql)

pers$pid <- pers$pid - 1
pers$home_id <- pers$home_id - 1
pers$day_id <- pers$day_id - 1
pers$day_id[is.na(pers$day_id)] <- -1
pers$undlycond[is.na(pers$undlycond)] <- -1

message("Preview for pers")
print(head(pers))
fwrite(pers %>% select(-undlycond),
       paste0(out_folder, "population-", cnt, ".txt"), sep = " ", scipen = 10)
fwrite(pers %>% select(-home_id, -day_id),
       paste0(out_folder, "comorbidity-", cnt, ".txt"), sep = " ", scipen = 10)

rm(hh_nonnh, loc, pers, wp)

## Household network (pretty much just take from the sql db)
sql <- "SELECT * FROM hh_network"
hh_network <- dbGetQuery(con, sql)
hh_network$locid1 <- hh_network$locid1 - 1
hh_network$locid2 <- hh_network$locid2 - 1

message("Preview for hh_network")
print(head(hh_network))
fwrite(hh_network, paste0(out_folder, "network-", cnt, ".txt"),
       sep = " ", col.names = F, scipen = 10)

rm(hh_network)

## Neighbourhood network (pretty much just take from the sql db)
# sql <- "SELECT * FROM nb_network"
# nb_network <- dbGetQuery(con, sql)
# nb_network$locid1 <- nb_network$locid1 - 1
# nb_network$locid2 <- nb_network$locid2 - 1
#
# message("Preview for nb_network")
# print(head(nb_network))
# fwrite(nb_network, paste0(out_folder, "network_nb-", cnt, ".txt"),
#        sep = " ", col.names = F, scipen = 10)
#
# rm(nb_network)

## Extracurricular activities (pretty much just take from the sql db)
if (extracurricular_mode) {
  sql <- "SELECT * FROM extracurr"
  ec <- dbGetQuery(con, sql)
  ec <- ec - 1

  message("Preview for public-activity")
  print(head(ec))
  fwrite(ec, paste0(out_folder, "public-activity-", cnt, ".txt"), sep = " ")

  ec_vec <- ec[,-1] %>%
    as.matrix() %>%
    as.vector()
  ec_vec <- unique(ec_vec) %>% sort()
  ec_df <- data.frame(ec = ec_vec)
  message("Preview for public-locations")
  print(head(ec_df))
  fwrite(ec_df, paste0(out_folder, "public-locations-", cnt, ".txt"), sep = " ",
         col.names = F)
}


#### Finish and kill the unpacked files
message("Removing unpacked files...")
dbDisconnect(con)
unlink(folder, recursive = T)
message("Done.")

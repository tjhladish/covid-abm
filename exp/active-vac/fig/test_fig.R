
.pkgs <- c("data.table", "ggplot2", "scales", "cabputils", "RSQLite")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

.args <- commandArgs(args = c(
  "~/Downloads/vpdr", # TODO some digest of this output
  file.path("covid-active-v7.sqlite"),
  file.path("fig", "vis_support.rda"),
  file.path("fig", "output", "inf_risk.png")
))

dt <- .args[1] |> list.files(pattern = "vpdr.*csv$", full.names = TRUE) |>
  setNames(nm = _) |>
  lapply(fread, drop = 1) |> rbindlist(idcol = "file") |>
  melt.data.table(id.vars = c("day", "file")) |>
  DT(, serial := gsub("^.*_(\\d+)\\.csv$", "\\1", file) |> as.integer())

dt$file <- NULL

conn <- dbConnect(SQLite(), .args[2])
serialkey <- conn |> dbGetQuery(
  "SELECT quar, pas_vac, act_vac, pas_alloc, act_alloc, inf_con, ppb_ctrl, min(serial), max(serial), count(*) FROM par WHERE season != 0.0 GROUP BY quar, pas_vac, act_vac, pas_alloc, act_alloc, inf_con, ppb_ctrl;"
) |> as.data.table()
dbDisconnect(conn)

.args[3] |> load()

translator <- dt[, serialkey[between(serial, `min(serial)`, `max(serial)`), .(quar, pas_vac, act_vac, pas_alloc, act_alloc, inf_con)], by=serial]

plt.dt <- dt[translator, on=.(serial), .(
  date = day + as.Date("2020-02-10"),
  variable, value,
     vtype = factor(c("none", "ring", "risk", "age")[act_vac+1]),
     alloc = factor(
       c("none", "LS", "MS", "HS", "USA")[ifelse(pas_vac == 0, act_alloc, pas_alloc)+1],
       levels = c("none", "LS", "MS", "HS", "USA"), ordered = TRUE
     ),
     quar,
     inf_con
   )
][inf_con != 2]

plt.dt[variable == "inf_risk", variable := "infrisk_tot"]
plt.dt[variable %like% "inf_risk", variable := gsub("inf_risk", "infrisk", variable) ]

plt.dt <- plt.dt[
  variable %like% "infrisk"
][, var := gsub("infrisk_", "", variable) |> factor(levels = c("nat", "vax", "tot"), ordered = TRUE) ]

p2 <- ggplot(plt.dt) +
  aes(x = date, y = value, color = vtype, linetype = factor(c("nonpi", "wquar")[quar+1])) +
  facet_grid(
    var ~ alloc, switch = "y",
    labeller = labeller(var = c(
      nat = "Prior Infection",
      vax = "Vaccination",
      tot = "All Sources"
    ))
  ) +
  geom_month_background(
    plt.dt[!is.nan(value)], by = c(row="var", col="alloc"), font.size = 3
  ) +
  geom_line() +
  theme_minimal() + theme(
    legend.position = "bottom", strip.placement = "outside"
  ) +
  scale_y_continuous(name = "Mean Infection Risk Given Exposure,\nReduced Due to ...") +
  scale_x_null() +
  scale_color_strategy(
    guide = guide_legend(title.position = "top", title.hjust = 0.5)
  ) + scale_linetype_quar()

store(p2, .args, width = 10, height = 6, bg = "white")

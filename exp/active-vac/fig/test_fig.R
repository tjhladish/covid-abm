
.pkgs <- c("data.table", "ggplot2", "cabputils")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

.args <- commandArgs(args = c(
  "~/Downloads", # TODO some digest of this output
  file.path("~", "Downloads", "serial_lookup.csv"), # TODO extract
  file.path("fig", "vis_support.rda"),
  "something.png"
))

dt <- .args[1] |> list.files(pattern = "vpdr.*csv$", full.names = TRUE) |>
  lapply(fread) |> rbindlist() |> melt.data.table(id.vars = c("serial", "day"))

serialkey <- .args[2] |> fread()
.args[3] |> load()

translator <- dt[, serialkey[between(serial, `min(serial)`, `max(serial)`), .(quar, pas_vac, act_vac, pas_alloc, act_alloc, inf_con)], by=serial]

dt[translator, on=.(serial), vtype := c("none", "ring", "risk", "age")[act_vac+1] ]

dt[variable == "inf_risk", variable := "infrisk_tot"]
dt[variable %like% "inf_risk", variable := gsub("inf_risk", "infrisk", variable) ]


plt.dt <- dt[
  variable %like% "infrisk"
][, var := gsub("infrisk_", "", variable) |> factor(levels = c("nat", "vax", "tot"), ordered = TRUE) ]

p2 <- ggplot(plt.dt) +
  aes(x = day + as.Date("2020-02-10"), y = value, color = vtype) +
  facet_grid(
    var ~ ., switch = "y",
    labeller = labeller(var = c(
      nat = "Infection Derived",
      vax = "Vaccine Derived",
      tot = "All Sources"
    ))
  ) +
  geom_month_background() +
  geom_line() +
  theme_minimal() + theme(
    legend.position = "bottom", strip.placement = "outside"
  ) +
  scale_y_continuous(

  ) +
  scale_x_null() +
  scale_color_strategy(
    guide = guide_legend(title.position = "top", title.hjust = 0.5)
  )

store(.args, p, width = 10, height = 6, bg = "white")


.pkgs <- c("data.table", "ggplot2", "cabputils")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

.args <- commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", "digest.rds"),
  file.path("fig", "process", "digest-key.rds"),
  file.path("fig", "output", "summary.png")
))

load(.args[1])

overdates <- as.Date(c("2021-05-31", "2021-11-30", "2022-03-31"))

scn.dt <- readRDS(.args[3])[inf_con == FALSE][, .(
  scenario, quar, alloc = fifelse(pas_vac, pas_alloc, act_alloc),
  act_vac
)]

dt <- readRDS(.args[2])[
  (outcome == "deaths") & (scenario %in% scn.dt$scenario)
]

q.dt <- dt |>
  DT(,c("c.value", "c.averted") := .(cumsum(value), cumsum(averted)), by=scenario) |>
  DT(date %in% overdates) |>
  quantile(j=.(c.value, c.averted, c.effectiveness), sampleby = "realization")


plt.dt <- q.dt[scn.dt, on = .(scenario), nomatch = 0]
plt.dt[date == overdates[1], variant := "alpha"]
plt.dt[date == overdates[2], variant := "delta"]
plt.dt[date == overdates[3], variant := "omicron"]

plt.dt$act_vac <- factor(plt.dt$act_vac, levels = c("ring", "none", "age", "risk"))

p <- ggplot(plt.dt) + aes(x=variant, color = act_vac, shape = quar) +
  geom_point(aes(y=qmed), data = \(dt) dt[quar == FALSE], position = position_dodge(width = 0.5)) +
  geom_point(aes(y=qmed), data = \(dt) dt[quar == TRUE], position = position_dodge(width = 0.5)) +
  facet_grid(measure ~ alloc, scales = "free_y", switch = "y") +
  theme_minimal() +
  theme(strip.placement = "outside") +
  scale_color_strategy(breaks = c("ring", "none", "age", "risk")) +
  scale_x_discrete("Post Variant Wave Outcome") +
  scale_y_continuous("Relative to Deaths ...")

store(obj = p, args = .args, width = 9, height = 5, bg = "white")

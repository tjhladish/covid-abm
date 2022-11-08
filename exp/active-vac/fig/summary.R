
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
  DT(,c("c.value", "c.averted") := .(cumsum(value), cumsum(averted)), by=.(scenario, realization)) |>
  DT(date %in% overdates) |>
  quantile(j=.(c.value, c.averted, c.effectiveness), sampleby = "realization", probs = qprobs(c(`90`=0.9)))


plt.dt <- q.dt[scn.dt, on = .(scenario), nomatch = 0]
plt.dt[date == overdates[1], variant := "alpha"]
plt.dt[date == overdates[2], variant := "delta"]
plt.dt[date == overdates[3], variant := "omicron"]

plt.dt$act_vac <- factor(plt.dt$act_vac, levels = c("ring", "none", "age", "risk"))

scale_shape_quar <- rejig(
  scale_shape_manual,
  name = "Extra NPI", labels = c(nonpi="None", wquar = "Quarantine Contacts"),
  values = c(wquar=17, nonpi=16),
  guide = guide_legend(title.position = "top", title.hjust = 0.5, order = 1)
)

p <- ggplot(plt.dt) + aes(
  x=variant, color = act_vac,
  shape = c("nonpi","wquar")[quar+1],
  linetype = c("nonpi","wquar")[quar+1]
) +
  geom_pointrange(
    aes(y=qmed, ymin=q90l, ymax=q90h),
    data = \(dt) dt[quar == FALSE],
    position = position_dodge(width = 0.5),
    size = .25
  ) +
  geom_pointrange(
    aes(y=qmed, ymin=q90l, ymax=q90h),
    data = \(dt) dt[quar == TRUE],
    position = position_dodge(width = 0.5),
    size = .25
  ) +
  facet_grid(
    measure ~ alloc, scales = "free_y", switch = "y",
    labeller = labeller(measure = c(
      c.value="Per 10K, Cumulative\nIncidence",
      c.averted="Per 10K, Cumulative\nAverted Incidence",
      c.effectiveness="Cumulative Relative\nAverted Incidence"
    ))
  ) +
  coord_cartesian(ylim = c(0, NA), expand = FALSE, clip = "off") +
  theme_minimal() +
  theme(
    strip.placement = "outside", legend.position = "bottom",
    panel.spacing.y = unit(1.5, "line"), panel.spacing.x = unit(1, "line"),
    legend.text = element_text(size = rel(.75))
  ) +
  scale_linetype_quar(
    guide = guide_legend(title.position = "top", title.hjust = 0.5, order = 1)
  ) + scale_shape_quar() +
  scale_color_strategy(breaks = c("ring", "none", "age", "risk")) +
  scale_x_discrete("Post Variant Wave Outcome") +
  scale_y_continuous("Relative to Deaths ...")

store(obj = p, args = .args, width = 9, height = 6, bg = "white")

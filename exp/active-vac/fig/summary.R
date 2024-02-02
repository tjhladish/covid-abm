
.pkgs <- c("data.table", "ggplot2", "cabputils")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

.args <- commandArgs(args = c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", "alt_eff.rds"),
  file.path("fig", "process", "digest-key.rds"),
  file.path("fig", "output", "summary.png")
))

load(.args[1])

overdates <- as.Date(c("2021-05-27", "2021-11-26", "2022-03-07"))

scn.dt <- readRDS(.args[3])[inf_con == FALSE][, .(
  scenario, quar, alloc = fifelse(pas_vac, pas_alloc, act_alloc),
  act_vac
)]

dt <- readRDS(.args[2])[
  (outcome %in% c("inf", "deaths")) & (scenario %in% scn.dt$scenario)
]

q.dt <- dt |>
  DT(,c("c.value", "c.averted") := .(cumsum(value), cumsum(averted)), by=setdiff(key(dt), "date")) |>
  DT(date %in% overdates) |>
  quantile(j=.(c.value, c.averted, c.effectiveness), sampleby = "realization", probs = qprobs(c(`90`=0.9)))


plt.dt <- q.dt[scn.dt, on = .(scenario), nomatch = 0]
ref.dt <- CJ(
  outcome = c("inf", "deaths"),
  date = overdates,
  measure = "c.effectiveness",
  qmed = 0, quar = FALSE, alloc = c("LS", "MS", "HS", "USA"),
  act_vac = "none"
)
plt.dt <- rbind(plt.dt, ref.dt, fill = TRUE)
plt.dt[date == overdates[1], variant := "alpha"]
plt.dt[date == overdates[2], variant := "delta"]
plt.dt[date == overdates[3], variant := "omicron"]

plt.dt$act_vac <- factor(plt.dt$act_vac, levels = c("ring", "none", "age", "risk"))

scale_shape_quar <- rejig(
  scale_shape_manual,
  name = "Extra NPI", labels = c(nonpi="None", wquar = "Quarantine Contacts"),
  values = c(wquar=21, nonpi=19),
  guide = guide_legend(title.position = "top", title.hjust = 0.5, order = 1)
)

p <- ggplot(plt.dt[measure == "c.effectiveness"]) + aes(
  x=variant, color = act_vac,
  shape = c("nonpi","wquar")[quar+1]
) +
  geom_pointrange(
    aes(y=qmed, ymin=q90l, ymax=q90h),
    data = \(dt) dt[quar == FALSE],
    position = position_dodge(width = 0.5),
    size = 0.5, stroke = 0
  ) +
  geom_pointrange(
    aes(y=qmed, ymin=q90l, ymax=q90h),
    data = \(dt) dt[quar == TRUE],
    position = position_dodge(width = 0.5), fill = "white",
    size = 0.4, stroke = 0.4
  ) +
  facet_grid(
    outcome ~ alloc, scales = "free_y", switch = "y",
    labeller = labeller(outcome = c(inf = "Infections", deaths = "Deaths"))
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    strip.placement = "outside", legend.position = "bottom",
    #panel.spacing.y = unit(1.5, "line"), panel.spacing.x = unit(1, "line"),
    legend.text = element_text(size = rel(.75)),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey")
  ) +
  # scale_linetype_quar(
  #   guide = guide_legend(title.position = "top", title.hjust = 0.5, order = 1)
  # ) +
  scale_shape_quar() +
  scale_color_strategy(
    breaks = c("ring", "none", "age", "risk"),
    values = c(none = "black", ring = "#fb6502", risk = "#3b90db", age = "#209033")
  ) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "Cumulative Effectiveness After Each Variant")

g <- ggplotGrob(p)
id <- which(g$layout$name == "guide-box")
g$layout[id, c("l","r")] <- c(1, ncol(g))
grDevices::png(tail(.args,1), width = 9, height = 6, bg = "white", units = "in", res = 300)
grid::grid.draw(g)
grDevices::dev.off()

#store(obj = g, args = .args, width = 9, height = 6, bg = "white")

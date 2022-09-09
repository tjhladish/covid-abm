
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("digest-doses.rds", "digest-key.rds", "alt_eff.rds")),
  file.path(c("covax_doses_COVAX_only.txt", "covax_doses_MIC_only.txt")),
  file.path("fig", "output", "cum_doses.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

#incovax <- fread(.args[5])

tar <- "doses"

scn.dt <- readRDS(.args[3])
addscn <- scn.dt[act_alloc == "none" & !(pas_alloc %in% c("none", "FL"))]
scns <- c(
  readRDS(.args[4])[, unique(scenario)],
  scn.dt[act_alloc == "none" & !(pas_alloc %in% c("none", "FL")), unique(scenario)]
)

#' comes key'd
doses.dt <- readRDS(.args[2])[
  eval(datefilter)
][scenario %in% scns, .(
  scenario, realization, date, c.value
)]

intscns <- doses.dt[, unique(scenario)]

plt.dt <- setkeyv(
  doses.dt[scn.dt[scenario %in% intscns], on=.(scenario)],
  union(key(doses.dt), colnames(scn.dt))
)

rm(doses.dt)
gc()

plt.dt[, talloc := factor(
  fifelse(
    pas_alloc == "none", as.character(act_alloc), fifelse(
    pas_alloc %in% c("FL", "FL+ring"), "FL+", as.character(pas_alloc)
  )), levels = c("COVAX", "MIC", "FL+"), ordered = TRUE
)][, qfac := factor(c("No Additional NPI", "Quarantine Contacts")[quar+1]) ]

p <- ggplot(plt.dt) + aes(
  x = date, y = c.value,
  color = act_vac,
  linetype = factor(c("unconditional", "conditional")[inf_con+1]),
  sample = realization
) +
  facet_nested(rows = vars(qfac), cols = vars(talloc)) +
  geom_month_background(
    plt.dt, by = c("qfac", "talloc"),
    font.size = 3, value.col = "c.value"
    # ,
    # ymax = plt.dt[, max(c.value), by=.(talloc)],
    # ymin = plt.dt[, min(c.value), by=.(talloc)]
  ) +
  stat_spaghetti(
    aes(alpha = after_stat(sampleN^-1))
  ) +
#  geom_hline(aes(yintercept=0, color = "none")) +
  scale_y_continuous(
    name = sprintf(
      "Cumulative Vaccination"
    )
  ) +
  scale_x_null() +
  scale_color_discrete(
    "Vaccine Program",
    breaks = c("risk", "ring", "none"),
    labels = c(
      ring="Infection-risk Prioritization",
      risk="Disease-risk Prioritization",
      none="Standard Only"
    )) +
  scale_linetype_manual("Vaccinate...", labels = c(conditional="Condtional on\nCase History", unconditional = "Unconditionally"), values = c(conditional="dashed", unconditional="solid")) +
  scale_alpha(range = c(0.01, 1)) +
  theme_minimal() +
  theme(
    legend.position = c(0, 1), legend.justification = c(0, 1), strip.placement = "outside"
  )

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")

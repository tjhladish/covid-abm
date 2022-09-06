
.pkgs <- c("data.table", "ggplot2", "scales", "cabputils")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("digest.rds", "digest-key.rds", "digest-ref.rds")),
  file.path("fig", "output", "averted.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

#' comes key'd
eff.dt <- readRDS(.args[2])[
  date > "2020-12-01"
][, .(
  scenario, realization, date, outcome, averted
)]

intscns <- eff.dt[, unique(scenario)]

scn.dt <- readRDS(.args[3])

plt.dt <- setkeyv(
  eff.dt[scn.dt[scenario %in% intscns], on=.(scenario)],
  union(key(eff.dt), colnames(scn.dt))
)

rm(eff.dt)
gc()

plt.dt[, scn := factor(sprintf(
  "%s%s\n%s\n%s",
  fifelse(quar,"q+, ",""),
  fifelse(inf_con,"[ic] ",""),
  fifelse(pas_vac, paste0("p=",pas_alloc,", "), ""),
  fifelse(act_vac != "none", paste0("a=",act_vac,":",act_alloc), "")
))]

asinh_trans <- function() {
  trans_new("asinh", asinh, inverse = sinh)
}

p <- ggplot(plt.dt) + aes(
  x = date, y = averted,
  color = outcome, sample = realization
) +
  facet_wrap(~scn) +
  geom_month_background(
    plt.dt, by = "scn",
    font.size = 3, value.col = "averted",
    ymax = plt.dt[, max(averted)],
    ymin = plt.dt[, min(averted)],
    ytrans = "asinh"
  ) +
  stat_spaghetti() +
  scale_y_averted(
    trans = "asinh",
    breaks = c(-200, -50, -10, -2, 0, 2, 10, 50, 200),
    limits = c(-200, 200)
  ) +
  scale_x_null() +
  theme_minimal() +
  theme(
    legend.position = "bottom", strip.placement = "outside"
  )

ggsave(tail(.args, 1), p, height = 10, width = 8, bg = "white")

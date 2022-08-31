
.pkgs <- c("data.table", "ggplot2", "scales", "ggh4x", "cabputils")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("digest.rds", "digest-key.rds", "digest-ref.rds")),
  file.path("fig", "output", "averted_inf.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

tar <- match.arg(
  tail(.args, 1) |> basename() |> gsub(pattern = "^.+_(.+)\\..+$", replacement = "\\1"),
  c("inf", "sev", "deaths")
)

#' comes key'd
eff.dt <- readRDS(.args[2])[
  date > "2020-12-01"
][
  outcome == tar
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

# plt.dt[, scn := factor(sprintf(
#   "%s%s\n%s\n%s",
#   fifelse(quar,"q+, ",""),
#   fifelse(inf_con,"[ic] ",""),
#   fifelse(pas_vac, paste0("p=",pas_alloc,", "), ""),
#   fifelse(act_vac != "none", paste0("a=",act_vac,":",act_alloc), "")
# ))]

asinh_trans <- function() {
  trans_new("asinh", asinh, inverse = sinh)
}

plt.dt[, pvac := factor(
  fifelse(pas_alloc %in% c("none", "COVAX", "MIC"), "Limited Passive", "High Passive")
)][, qfac := factor(c("No Q", "Q")[quar+1]) ]

p <- ggplot(plt.dt) + aes(
  x = date, y = averted,
  color = interaction(act_vac, act_alloc, drop = TRUE),
  linetype = factor(c("unconditional", "conditional")[inf_con+1]),
  sample = realization
) +
  facet_nested(rows = vars(qfac), cols = vars(pvac, pas_alloc)) +
  geom_month_background(
    plt.dt, by = c("qfac", "pvac", "pas_alloc"),
    font.size = 3, value.col = "averted",
    ymax = plt.dt[, max(averted)],
    ymin = plt.dt[, min(averted)],
    ytrans = "asinh"
  ) +
  stat_spaghetti(
    aes(alpha = after_stat(sampleN^-1))
  ) +
  scale_y_averted(
    name = sprintf(
      "Per 10k, Averted Incidence\nof %s",
      switch(tar, inf = "Infection", sev = "Severe Disease", deaths = "Death", stop())
    ),
    trans = "asinh"
  ) +
  scale_x_null() +
  scale_color_discrete("Active Vax.") +
  scale_linetype_discrete("Conditional Vax.") +
  scale_alpha(range = c(0.01, 1)) +
  theme_minimal() +
  theme(
    legend.position = "bottom", strip.placement = "outside"
  )

ggsave(tail(.args, 1), p, height = 6, width = 10, bg = "white")

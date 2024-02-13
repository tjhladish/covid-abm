
.pkgs <- c("data.table", "ggplot2", "ggh4x", "scales", "cabputils")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c("digest-doses.rds", "digest-key.rds")),
  file.path("fig", "output", "doses.png")
) else commandArgs(trailingOnly = TRUE)

load(.args[1])

#' comes key'd
#' TODO: identify proper reference point?
doses.dt <- readRDS(.args[2])[eval(datefilter)]

scn.dt <- readRDS(.args[3])[inf_con == FALSE][eval(seasfilter)]

plt.dt <- setkeyv(
  doses.dt[scn.dt, on=.(scenario), nomatch = 0],
  union(key(doses.dt), colnames(scn.dt))
)[eval(seasfilter)]

plt.qs <- quantile(
  plt.dt,
  j = .(c.value), sampleby = "realization",
  probs = qprobs(c(`90`=.9), mid = TRUE, extent = FALSE)
)[, talloc := factor(
  fifelse(
    pas_vac == FALSE,
    as.character(act_alloc), as.character(pas_alloc)
  ), levels = c("LS", "MS", "HS", "USA"), ordered = TRUE
)][, qfac := factor(c("No Additional NPI", "Quarantine Contacts")[quar+1]) ][,
  outcome := "vaccine"
]

p <- allplot(
  plt.qs, yl = NULL, withRef = FALSE
) + theme(strip.text.y = element_text(size = rel(1.25))) +
  geom_line(aes(y=qmed), data = \(dt) subset(dt, pas_vac == TRUE & quar == FALSE), show.legend = FALSE, size = 0.25)

store(p, .args, height = 3, width = 10, bg = "white")

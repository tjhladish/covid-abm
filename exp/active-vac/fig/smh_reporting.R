
.pkgs <- c("data.table", "scales", "ggplot2", "ggrepel", "patchwork", "cabputils", "geomtextpath")
.pkgs |> sapply(require, character.only = TRUE) |> all() |> stopifnot()

.args = commandArgs(
  trailingOnly = TRUE, c(
  file.path("fig", "vis_support.rda"),
  file.path("fig", "process", c(
    "validation.rds",
    "detection.rds"
  )),
  file.path("fig", "smh_reporting.png")
))

intfilter <- if (interactive()) expression(realization < 10) else expression(realization >= 0)

load(.args[1])
d <- readRDS(.args[2])[date <= vendday][eval(intfilter)]
ref.day0 <- d[, min(date)]
detect.dt <- readRDS(.args[3])[, date := day + ref.day0][date <= vendday]

conserved <- list(
  scale_x_null(),
  scale_shape_measure(),
  theme_minimal(),
  theme(text = element_text(face = "bold")),
  scale_alpha_continuous(guide = "none", range = c(0.05, 1))
)

geom_grid <- function(
  ys, FUN = as.character
) {
  geom_texthline(
    aes(yintercept = ys, label = label),
    data = data.table(ys = ys, label = FUN(ys)), inherit.aes = FALSE,
    gap = TRUE, hjust = 0.325, color = "grey75", fontface = "bold"
  )
}

geom_month_end <- function(
    data,
    col.cycle = c(off = NA, on = alpha("grey95", 0.75)),
    m.labels = m.abb,
    font.size = 8, font.face = "bold",
    ytrans = "identity",
    datafn = tsref.dt,
    m.y = 0.0, y.y = 0.3,
    text.col = rep("darkgrey", length(col.cycle)),
    ...
) {
  names(text.col) <- c(tail(names(col.cycle), -1), names(col.cycle)[1])
  text.col[is.na(text.col)] <- "white"
  dt <- datafn(data, ...)
  dt[, fill := col.cycle[mtype] ]
  dt[, col := text.col[mtype] ]

  xformer <- scales::as.trans(ytrans)

  xform <- function(ymax, ymin, dropto) {
    xformer$inverse(
      (xformer$transform(ymax)-xformer$transform(ymin))*dropto + xformer$transform(ymin)
    )
  }

  # ggplot isn't great on replicating these across facets
  # would be preferrable to use `annotate`, and let backend recycle fills, but
  # it won't, so have to tell this how many facets there will be
  list(
    geom_rect(
      mapping = aes(xmin = start - 0.5, xmax = end + 0.5, ymin = if (xformer$name %like% "log") 0 else -Inf, ymax = if (xformer$name %like% "prob") 1 else Inf),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, fill = dt$fill
    ),
    geom_text(
      mapping = aes(x = mid, y = xform(ymax, ymin, m.y), label = m.abb[mon]),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, color = dt$col,
      size = font.size,
      vjust = if (y.y > 0) "bottom" else "top",
      fontface = font.face
    ),
    geom_text(
      mapping = aes(x = start+5, y = xform(ymax, ymin, y.y), label = yr),
      data = dt[yshow == TRUE],
      inherit.aes = FALSE, show.legend = FALSE, color = dt[yshow == TRUE]$col,
      size = font.size, hjust = 0, fontface = font.face,
      vjust = if (y.y > 0) "top" else "bottom"
    )
  )
}

p.bot <- ggplot() +
  geom_month_end(
    d[, .(realization, date, value = cov1) ],
    ymax = 1, ymin = 0, m.y = -0.05, y.y = -0.6,
    font.size = 6, text.col = rep("black", 2)
  ) + conserved + theme(
    axis.text = element_blank(), axis.title = element_blank(),
    panel.grid = element_blank(), plot.margin = margin()
  ) + coord_cartesian(ylim=c(-0.6,0))

geom_month_bars <- function(
    data,
    col.cycle = c(off = NA, on = alpha("grey95", 0.75)),
    ytrans = "identity",
    datafn = tsref.dt,
    ...
) {
  dt <- datafn(data, ...)
  dt[, fill := col.cycle[mtype] ]

  xformer <- scales::as.trans(ytrans)

  xform <- function(ymax, ymin, dropto) {
    xformer$inverse(
      (xformer$transform(ymax)-xformer$transform(ymin))*dropto + xformer$transform(ymin)
    )
  }

  # ggplot isn't great on replicating these across facets
  # would be preferrable to use `annotate`, and let backend recycle fills, but
  # it won't, so have to tell this how many facets there will be
  list(
    geom_rect(
      mapping = aes(xmin = start - 0.5, xmax = end + 0.5, ymin = if (xformer$name %like% "log") 0 else -Inf, ymax = if (xformer$name %like% "prob") 1 else Inf),
      data = dt, inherit.aes = FALSE, show.legend = FALSE, fill = dt$fill
    )
  )
}

geom_liner <- function(datafn) geom_text(
  aes(label = lbl, vjust = vj, hjust = hj, shape = NULL),
  data = datafn,
  fontface = "bold", size = 5
)

p.core <- function(
  dt, ymin = NA, ymax = NA, ytrans = "identity", gridy
) ggplot(dt) +
  aes(date, value) +
#  geom_grid(gridy) +
  geom_month_bars(dt, ymin = ymin, ymax = ymax, ytrans = ytrans) +
  theme_minimal() + theme(text = element_text(face = "bold")) +
  scale_x_null()

det.dt <- prepare(detect.dt[, .(
  realization = 1, date = day + ref.day0, asymp, mild, severe, crit
)])

setnames(det.dt, "measure", "outcome")
det.dt[, measure := "detection"]

p.detect <- p.core(
  det.dt, ymin = 0, ymax = 1, gridy = c(0.25, 0.5, 0.75)
) + aes(linetype = outcome) +
  geom_line() +
  scale_y_fraction(name = "Probability of Detection") +
  geom_text(mapping = aes(label=lab, linetype = NULL, shape = NULL, hjust = align), data = data.table(
    lab = c("Asymptomatic", "Mild", "Severe", "Critical"),
    date = as.Date(c("2020-07-15", "2020-05-10", "2020-04-15", "2020-08-22")),
    value = c(.18, .22, .65, .9),
    align = c(0.5, 1, 0.5, 0.5),
    measure = "detection"
  ), size = 6) +
  scale_linetype_manual(
    name = NULL, breaks = rev(c("asymp", "mild", "severe", "crit")),
    labels = c(asymp = "Asymptomic", mild = "Mild", severe = "Severe", crit = "Critical"),
    values = c(asymp = "dotted", mild = "dotdash", severe = "longdash", crit = "solid")
  ) +
  theme(
    legend.position = "none", panel.grid.major.y = element_blank(),
    plot.margin = margin(t=0.5, b=0, unit = "line"),
    axis.text.y = element_text(face="bold", color = "black", size = 12),
    axis.title.y = element_text(face="bold", color = "black", size = 12)
  ) + coord_cartesian(clip = "off")

p.res <- p.detect + p.bot +
  plot_layout(ncol = 1, heights = c(1, 0.15)) +
  theme(
    plot.tag.position = c(0, 0.97),
    axis.title = element_blank(), axis.text = element_blank(),
    plot.tag = element_text(size = rel(1.65), hjust = 0, vjust = 1)
  )

store(p.res, .args, width = 14, height = 4, bg = "white")

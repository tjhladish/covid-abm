
require(data.table)
require(ggplot2)

.args <- if (interactive()) c("test.png") else commandArgs(trailingOnly = TRUE)

ve.res <- 0.05
ve.step <- 0.2

# TODO move contouring logic in here?
dt <- data.table(expand.grid(R=seq(1,3,length.out = 101), VE=seq(.9, .1, by=-ve.res)))
dt[, lower := (1-(1/R))/VE ]
dt[, upper := (1-(1/R))/(VE-ve.res) ]
dt[VE == .1, upper := 1]

legend.inset <- 0.05
ve.contours <- seq(.1, .9, by=ve.step)

#' TODO probably superior to annotate the contour lines
#'  and provide a color bar instead of a legend
p <- ggplot(dt) + aes(R, group = VE) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = VE), alpha = 0.4) +
  geom_line(
    aes(y=lower, color = VE),
    data = function(dt) dt[abs(round((VE - 0.1)/ve.step)-(VE - 0.1)/ve.step) < ve.res/2]
  ) +
  theme_minimal() + theme(
    legend.position = c(1-legend.inset, 0+legend.inset),
    legend.justification = c(1, 0)
  ) +
  coord_cartesian(ylim=c(0, 1), expand = FALSE) +
  scale_color_continuous(
    expression(italic(VE)),
    breaks = ve.contours,
    labels = scales::label_percent(),
    guide = "legend", aesthetics = c("color","fill")
  ) +
  scale_x_continuous(expression("Effective reproductive number, "*italic(R)[eff])) +
  scale_y_continuous(
    expression("Critical vaccination coverage, "*italic(V)[C]),
    breaks = seq(0,1,by=0.2),
    labels = scales::label_percent()
  )

ggsave(tail(.args, 1), p, width = 5, height = 5, units = "in", dpi = 300, bg = "white")

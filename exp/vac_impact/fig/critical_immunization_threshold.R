
require(data.table)
require(ggplot2)

.args <- if (interactive()) c("test.png") else commandArgs(trailingOnly = TRUE)

ve.res <- 0.01
ve.step <- 0.2

# TODO move contouring logic in here?
dt <- data.table(
  expand.grid(R=seq(1,3,length.out = 101), VE=seq(ve.res, 1-ve.res, by=ve.res))
)
dt[, c("lower","upper","susc") := {
  f <- (1-(1/R))
  .(f/VE, f/(VE-ve.res), "full")
}]

pR <- 0.3

plot.dt <- rbind(
  dt,
  copy(dt)[, c("lower","upper","susc") := {
    f <- (1-pR-(1/R))
    .(f/(VE - VE*pR), f/((VE-ve.res) - (VE-ve.res)*pR), "part")
  }]
)

plot.dt[VE == min(VE), upper := 1]

legend.inset <- 0.01
ve.contours <- seq(.1, .9, by=ve.step)

#' TODO probably superior to annotate the contour lines
#'  and provide a color bar instead of a legend
p <- ggplot(plot.dt) + aes(R, group = VE) +
  facet_grid(. ~ susc, labeller = labeller(susc = function(l) {
    fifelse(l == "full", "Fully Susceptible", sprintf("%i%% Recovered", pR*100))
  })) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = VE), alpha = 0.4) +
  geom_line(
    aes(y=lower, color = VE),
    data = function(dt) dt[abs(round((VE - 0.1)/ve.step)-(VE - 0.1)/ve.step) < ve.res/2]
  ) +
  theme_minimal() + theme(
    legend.position = c(1-legend.inset, 0+legend.inset),
    legend.justification = c(1, 0),
    legend.direction = "horizontal",
    panel.spacing.x = unit(2, "lines"),
    legend.key.width = unit(2, "lines")
  ) +
  coord_cartesian(ylim=c(0, 1), expand = FALSE) +
  scale_color_continuous(
    expression(italic(VE)),
    breaks = ve.contours,
    labels = scales::label_percent(),
    aesthetics = c("color","fill"),
    guide = guide_colorbar(title.position = "top")
  ) +
  scale_x_continuous(expression("Basic reproductive number, "*italic(R)[0])) +
  scale_y_continuous(
    expression("Critical vaccination coverage, "*italic(V)[C]),
    breaks = seq(0,1,by=0.2),
    labels = scales::label_percent()
  )

ggsave(tail(.args, 1), p, width = 5, height = 5, units = "in", dpi = 300, bg = "white")

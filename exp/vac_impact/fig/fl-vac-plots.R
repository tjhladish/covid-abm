
library(data.table)
library(ggplot2)
library(ggh4x)
library(grid)
library(scales)

.mutation <- FALSE
.args <- if (interactive()) c(
  "vaceff.rds",
  "rt.rds",
  "plottools.rda",
  as.integer(.mutation),
  sprintf("fl-vac-mutation-%i.png", as.integer(.mutation))
) else commandArgs(trailingOnly = TRUE)

muttar <- as.integer(tail(.args, 2)[1])
vac.impact <- readRDS(.args[1])[mutation == muttar]
rt.dt <- readRDS(.args[2])
load(.args[3])

#' TODO read from some config file?
totalpop_10k <- 375.5 # Adjust accordingly

int.start <- vac.impact[which.max(averted > 0), min(date)]

datebg <- c(month_background_scale(), event_vline(vline.date = int.start))

ind.alpha = 0.025

plot.mlt <- rbind(melt(
  vac.impact,
  id.vars = c("realization", "date", "var"),
  measure.vars = c("value", "caverted", "ceff")
)[, intervention := TRUE ],
  vac.impact[,.(realization, date, var, variable = "value", value = i.value, intervention = FALSE)]
)


inc.dt <- plot.mlt[variable == "value"][, value := value/totalpop_10k ]

plot.dt <- rbind(
  rt.dt,
  inc.dt
)

levels(plot.dt$variable) <- c(
  expression(R[t]),
  expression("Incidence"~"per"~"10k"),
  "Cum. Averted per 10k",
  "Cum. Effectiveness"
)

percap.plot <- ggplot(
  plot.dt[(intervention == FALSE) | date >= int.start]
) + aes(date, value, color = intervention) +
  facet_grid2(
    variable ~ var,
    scales = "free_y", independent = "y",
    switch = "both",
    labeller = labeller(
      variable = label_parsed,
      var = c(c="Symptomatic Infections", d="Deaths")
    ), drop = TRUE
  ) +
  datebg +
  geom_hline(
    aes(yintercept = 1),
    data = function(dt) dt[variable == "R[t]"],
    linetype = "dashed",
    alpha = 0.5
  ) +
  geom_line(aes(group = interaction(variable, realization, intervention, var)), alpha = ind.alpha) +
  geom_line(
    data = function(dt) dt[,
      .(value = median(value)),
      by=.(variable, var, intervention, date)
    ],
    lineend = "round"
  ) +
  coord_cartesian(expand = FALSE, clip = "off") +
  scale_y_continuous(name=NULL) +
  theme_minimal() + theme(
    strip.placement = "outside",
    legend.position = c(0.75, 0.75),
    legend.justification = c(0.5, 0.5),
    panel.grid.major.x = element_blank(),
    panel.spacing.y = unit(1, "lines")
  ) +
  scale_color_intervention()

epi.grob <- ggplotGrob(percap.plot)
idx <- which(epi.grob$layout$name %in% c("panel-2-1", "axis-l-1-2"))
for (i in idx) epi.grob$grobs[[i]] <- nullGrob()

ggsave(tail(.args, 1), epi.grob, bg = "white", height = 3, width = 6)

# eff.plot <- ggplot(
#   plot.mlt[variable == "c.eff"]
# ) + aes(date, value, color = intervention) +
#   facet_grid2(
#     variable ~ var,
#     switch = "y", scales = "free", independent = "y",
#     labeller = labeller(
#       var = c(c="Symptomatic Infections", d="Deaths"),
#       variable = c(
#         value="Incidence per 10k",
#         c.averted="Cum. Averted\nIncidence per 10k",
#         c.eff="Cum. Effectiveness"
#       )
#     )
#   ) +
#   datebg() +
#   geom_line(aes(group = interaction(realization, intervention)), alpha = ind.alpha) +
#   geom_line(
#     data = function(dt) dt[,
#                            .(value = median(value)),
#                            by=.(variable, var, intervention, date)
#     ]
#   ) +
#   scale_y_continuous(name=NULL) +
#   theme_minimal() + theme(
#     strip.placement = "outside",
#     strip.background.x = element_blank(),
#     strip.text.x = element_blank()
#   ) +
#   coord_cartesian(ylim = c(NA, 1), expand = FALSE) +
#   intervention_scale

# (Rt.plot + guide_area() + percap.plot + eff.plot) + plot_layout(
#   guides = "collect",
#   design="
# AB
# CC
# CC
# DD
# ")

# percap.plot + theme(
#   
# )
# 
# (Rt.plot + guide_area() + percap.plot) + plot_layout(
#   guides = "collect",
#   design="
# AB
# CC
# ")
# 
# plot_curves <- function (epc, param, legend = F, hline = NULL, ...) {
# 
#   flag <- ifelse(startsWith(param, "Eff") | param == "case_avert" |
#                    param == "death_avert",
#                  "n", "l")
# 
#   cond <- epc$vac == 0
#   plot(epc[cond, c("date", param)],
#        type = "n",
#        col = c("black"),
#        lwd = 2,
#        axes = F, ...)
#   axis(side = 1,
#        at = months,
#        labels = format(months, "%b-%y"),
#        cex = 0.75)
#   axis(side = 2)
# 
#   cond <- epc$vac == 1
#   lines(epc[cond, c("date", param)],
#         col = c("dodgerblue"),
#         lwd = 2)
# 
#   if (flag == "l") {
#     cond <- epc$vac == 0
#     lines(epc[cond, c("date", param)],
#           col = c("black"),
#           lwd = 2)
#   }
# 
#   if (legend) {
#     legend("right",
#            legend = c("No Vaccine", "Vaccine"),
#            col = c("black","dodgerblue"),
#            lty = 1,
#            lwd = 2,
#            bty = "n",
#            xpd = NA)
#   }
# 
#   if (!is.null(hline)) abline(h = hline, lty = 3, lwd = 2, xpd = F)
# }
# 
# #### Calculate Effectiveness, Averted and Stuff
# ## Read from sqlite
# sql <- paste0("SELECT vac, realization, mutation, serial FROM par")
# par <- dbGetQuery(con, sql)
# 
# ## Joining table
# dat <- par %>%
#   left_join(met, by = "serial")
# 
# #### Pivot longer ----
# ## Using vac == 0 as control for all other scenarios
# epi_long <- dat %>%
#   pivot_longer(-c(vac, mutation, realization, serial),
#                names_to = c(".value", "week"),
#                names_pattern = "([a-zA-Z]+)([0-9]+)")
# 
# epi_long <- epi_long %>%
#   mutate(incd_c = c / totalpop_10k,
#          incd_d = d / totalpop_10k,
#          week = as.numeric(week))
# 
# #### Calculate effectiveness and case averted ----
# ### Eff1 = Cumulative from beginning
# ### Eff2 = Cumulative since vaccination
# ### We can actually do realization matching here
# 
# ## Add incidence of control
# epi_long1 <- epi_long %>%
#   group_by(mutation, realization, week) %>%
#   mutate(incd_c_control = incd_c[vac == 0],
#          incd_d_control = incd_d[vac == 0]) %>%
#   arrange(realization, week)
# 
# ## For each trt combination and realization,
# ## calculate effectiveness & case averted
# start_wk <- 47
# epi_long1 <- epi_long1 %>%
#   group_by(vac, mutation, realization) %>%
#   mutate(Eff1_c = calc_eff(incd_c, incd_c_control, 1),
#          Eff2_c = calc_eff(incd_c, incd_c_control, start_wk),
#          Eff1_d = calc_eff(incd_d, incd_d_control, 1),
#          Eff2_d = calc_eff(incd_d, incd_d_control, start_wk),
#          case_avert = cumsum(incd_c_control) - cumsum(incd_c),
#          death_avert = cumsum(incd_d_control) - cumsum(incd_d))
# 
# ## Take median of all metrics along all realization
# epi_long1 <- epi_long1 %>%
#   group_by(vac, mutation, week) %>%
#   summarise(incd_c = median(incd_c),
#             incd_d = median(incd_d),
#             Rt = median(Rt),
#             Eff1_c = median(Eff1_c),
#             Eff2_c = median(Eff2_c),
#             Eff1_d = median(Eff1_d),
#             Eff2_d = median(Eff2_d),
#             case_avert = median(case_avert),
#             death_avert = median(death_avert))
# 
# epi_long1 <- epi_long1 %>%
#   left_join(weeks_df, by = "week")
# 
# 
# for (m in 0:1) {
#   epc <- epi_long1 %>%
#     ungroup() %>%
#     filter(mutation == m)
# 
#   outfile <- paste0("fig/fl-vac-mutation-", m, ".png")
# 
#   ## Plot start
#   png(outfile, width = 20, height = 16, units = "cm",
#       res = 200)
# 
#   par(mfrow = c(4, 2),
#       mar = c(2, 4, 2, 1),
#       oma = c(0, 0, 3, 0))
# 
#   plot_curves(epc, param = "Rt", legend = F, hline = 1,
#               xlab = "",
#               ylab = "",
#               # xlim = c(0, max_week - 1),
#               ylim = c(min(epi_long1$Rt), max(epi_long1$Rt)),
#               main = "Rt")
#   plot.new()
#   legend("center",
#          legend = c("No Vaccine", "Vaccine"),
#          col = c("black","dodgerblue"),
#          lty = 1,
#          lwd = 2,
#          bty = "n",
#          xpd = NA,
#          cex = 1.3)
#   plot_curves(epc, param = "incd_c", legend = F,
#               xlab = "",
#               ylab = "",
#               # xlim = c(0, max_week - 1),
#               ylim = c(0, max(epi_long1$incd_c)),
#               main = "Case incidence per 10k")
#   plot_curves(epc, param = "incd_d", legend = F,
#               xlab = "",
#               ylab = "",
#               # xlim = c(0, max_week - 1),
#               ylim = c(0, max(epi_long1$incd_d)),
#               main = "Death incidence per 10k")
#   plot_curves(epc, param = "Eff2_c", legend = F, hline = 0,
#               xlab = "",
#               ylab = "",
#               # xlim = c(0, max_week - 1),
#               #ylim = range(epi_long1$Eff2_c),
#               ylim = c(0,1),
#               main = "Cumu. effectiveness (case)")
#   plot_curves(epc, param = "Eff2_d", legend = F, hline = 0,
#               xlab = "",
#               ylab = "",
#               # xlim = c(0, max_week - 1),
#               #ylim = range(epi_long1$Eff2_d),
#               ylim = c(0,1),
#               main = "Cumu. effectiveness (death)")
#   plot_curves(epc, param = "case_avert", legend = F, hline = 0,
#               xlab = "",
#               ylab = "",
#               # xlim = c(0, max_week - 1),
#               ylim = range(epi_long1$case_avert),
#               main = "Cumu. case averted per 10k")
#   plot_curves(epc, param = "death_avert", legend = F, hline = 0,
#               xlab = "",
#               ylab = "",
#               # xlim = c(0, max_week - 1),
#               ylim = range(epi_long1$death_avert),
#               main = "Cumu. death averted per 10k")
# 
#   main <- ifelse(m == 0, "No mutant strain", "With mutant strain")
# 
#   # mtext("Week after first introduction", 1, 1, outer = T)
#   mtext(main, 3, 0, outer = T)
# 
#   dev.off()
# }
# 
# dbDisconnect(con)

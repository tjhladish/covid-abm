
.pkgs <- c("data.table", "ComplexUpset", "ggsci", "ggplot2", "patchwork", "plyr", "gridExtra", "grid")

stopifnot(all(sapply(.pkgs, require, character.only = TRUE)))

sim_pop = "pseudo-300K"
# sim_pop = "sim_pop-florida"

if (interactive()) { setwd("~/documents/work/covid-abm/exp/active-vac/fig") }

#' assumes R project at the experiment root level
.args <- if (interactive()) c(
  file.path("..", "..", "..", "pop", sim_pop, paste0("locations-", sim_pop, ".txt")),
  file.path("..", "..", "..", "pop", sim_pop, paste0("public-locations-", sim_pop, ".txt")),
  file.path("..", "..", "..", "pop", sim_pop, paste0("population-", sim_pop, ".txt")),
  file.path("..", "..", "..", "pop", sim_pop, paste0("network-", sim_pop, ".txt")),
  file.path("manuscript", "supp_locs_by_type.png"),
  file.path("manuscript", "supp_home_size_distr.png"),
  file.path("manuscript", "supp_work_size_distr.png"),
  file.path("manuscript", "supp_home_degree_distr.png")
) else commandArgs(trailingOnly = TRUE)

output_path = file.path('.', 'manuscript')
dir.create(output_path, recursive = TRUE)

locs.dt = fread(.args[1])[locid %in% fread(.args[2])[,V1], patronized := 1][is.na(patronized), patronized := 0]

home_size.dt = fread(.args[3])[, .N, by = .(home_id)]
work_size.dt = fread(.args[3])[, .N, by = .(day_id)]

network.dt = fread(.args[4])
edge1 = network.dt[, .(N1 = .N), by = .(locid = V1)]
edge2 = network.dt[, .(N2 = .N), by = .(locid = V2)]
edge.dt = edge1[edge2, on = .(locid)][is.na(N1), N1 := 0][is.na(N2), N2 := 0][, num_edges := N1 + N2]

# set.list = list(
#   `Essential` = locs.dt[essential == "y", locid],
#   `Non-essential` = locs.dt[essential == "n", locid],
#   `Patronizable` = locs.dt[transmission != "N", locid],
#   `Non-patronizable` = locs.dt[transmission == "N", locid],
#   `Normal risk` = locs.dt[transmission != "H", locid],
#   `High risk` = locs.dt[transmission == "H", locid]
# )
# 
# set.colors = c("#3B4992FF", "#EE0000FF", "#008B45FF", "#008B45FF", "#EE0000FF", "#3B4992FF")

upset.dt = locs.dt[type == 'w', .(
  `Essential` = ifelse(essential == "y", TRUE, FALSE),
  `Non-essential` = ifelse(essential == "n", TRUE, FALSE),
  `Patronizable` = ifelse(transmission != "N", TRUE, FALSE),
  `Non-patronizable` = ifelse(transmission == "N", TRUE, FALSE),
  `Normal risk` = ifelse(transmission != "H", TRUE, FALSE),
  `High risk` = ifelse(transmission == "H", TRUE, FALSE)
)]

set.cols = list(
  "Essential" = "#3d8bffFF",
  "Non-essential" = "#9ec5ffFF",
  "Patronizable" = "#f68080FF",
  "Non-patronizable" = "#EE0000FF",
  "Normal risk" = "#008B45FF",
  "High risk" = "#80c5a2FF"
)

# set.cols = list(
#   "Essential" = "#3d8bffFF",
#   "Non-essential" = "#3d8bffFF",
#   "Patronizable" = "#EE0000FF",
#   "Non-patronizable" = "#EE0000FF",
#   "Normal risk" = "#008B45FF",
#   "High risk" = "#008B45FF"
# )

p1 = upset(
  upset.dt,
  colnames(upset.dt),
  name = "Location type",
  width_ratio = 0.25,
  matrix=(
    intersection_matrix(geom=geom_point(shape='circle filled', size=3))
    + scale_color_manual(
      values=c(
        "Essential" = set.cols[["Essential"]],
        "Non-essential" = set.cols[["Non-essential"]],
        "Patronizable" = set.cols[["Patronizable"]],
        "Non-patronizable" = set.cols[["Non-patronizable"]],
        "Normal risk" = set.cols[["Normal risk"]],
        "High risk" = set.cols[["High risk"]]
      ),
      guide = "none"
    )
    + theme(legend.position = "none")
  ),
  queries=list(
    upset_query(set = "Essential", fill = set.cols[["Essential"]]),
    upset_query(set = "Non-essential", fill = set.cols[["Non-essential"]]),
    upset_query(set = "Patronizable", fill = set.cols[["Patronizable"]]),
    upset_query(set = "Non-patronizable", fill = set.cols[["Non-patronizable"]]),
    upset_query(set = "Normal risk", fill = set.cols[["Normal risk"]]),
    upset_query(set = "High risk", fill = set.cols[["High risk"]])
  )
)

# change desnity to TRUE
p2 = ggplot(home_size.dt[home_id %in% locs.dt[type == "h", locid]]) +
  aes(x = N) +
  geom_histogram(
    aes(y = ..count../sum(..count..)),
    binwidth = 1,
    fill = "lightgray",
    color = "black"
  ) +
  # scale_x_continuous(trans = "log10") +
  theme_light() +
  labs(x = "Household size", y = "")

# log transform the X axis
p3 = ggplot(work_size.dt[day_id %in% locs.dt[type == "w", locid]]) +
  aes(x = N) +
  geom_histogram(
    aes(y = ..count../sum(..count..)),
    binwidth = 0.25,
    fill = "lightgray",
    color = "black"
  ) +
  scale_x_continuous(trans = "log10") +
  theme_light() +
  labs(x = "Workplace size", y = "")

p4 = ggplot(edge.dt) +
  aes(x = num_edges) +
  geom_histogram(
    aes(y = ..count../sum(..count..)),
    binwidth = 1,
    fill = "lightgray",
    color = "black"
  ) +
  # scale_x_continuous(trans = "log10") +
  theme_light() +
  labs(x = "Inter-household degree", y = "")


png(filename = tail(.args, 4)[1], height = 2400, width = 4800, res = 480)
p1
dev.off()

ggsave(tail(.args, 4)[2], p2, height = 8, width = 10, bg = "white", dpi = 300)
ggsave(tail(.args, 4)[3], p3, height = 8, width = 10, bg = "white", dpi = 300)
ggsave(tail(.args, 4)[4], p4, height = 8, width = 10, bg = "white", dpi = 300)

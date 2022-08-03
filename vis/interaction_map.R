rm(list=ls())
setwd("~/work/covid-abm/pop/pseudo-300K")
library('data.table')
library('RColorBrewer')

p = fread('population-pseudo-300K.txt')
l = fread('locations-pseudo-300K.txt')
n = fread('network-pseudo-300K.txt')
a = fread('public-activity-pseudo-300K.txt')

net = n[l[,c('locid','x','y','compliance')], on = c(V1 = "locid"), nomatch=0]
net = net[l[,c('locid','x','y','compliance')], on = c(V2 = "locid"), nomatch=0]
names(net) = c('l1', 'l2', 'x1', 'y1', 'c1','x2', 'y2', 'c2')

net$type = 'social'
net$r = 1 - pmax(net$c1, net$c2)

houses = l[type=='h',]
houses$bin = findInterval(1 - houses$compliance, seq(0,1,0.01))
net$bin = findInterval(net$r, seq(0,1,length.out=256))

xlim = c(-82.2, -82.194)
ylim = c(29.0, 29.006)

plotter = function(base_filename, snet, lcols, as.png=T) {
  if(as.png) png(paste0(base_filename, '.png'), width=2000, height=2000, res=300)
  par(mar=c(3,3,1,1))
  plot(snet$x1, snet$y1, type='n', xlim=xlim, ylim=ylim, xlab='', ylab='')
  segments(snet$x1, snet$y1, snet$x2, snet$y2, col=lcols)
  pal = colorRampPalette(c("blue", "red"))
  points(houses$x, houses$y, pch=19, col=pal(100)[houses$bin])
  if(as.png) dev.off()
}

# is either end point of the edge within the plot limits?
plot_net = net[(x1 %between% xlim & y1 %between% ylim) | (x2 %between% xlim & y2 %between% ylim)]

subnet = plot_net[runif(nrow(plot_net)) < 0.5] #subset randomly to reduce density
line_cols = paste0('#000000', as.hexmode(subnet$bin)) # opacity = riskiness
plotter('interaction_map_shaded', subnet, line_cols)

# let's try to add more stuff
plotter('', subnet, line_cols, as.png=F)
points(l[type=='w',.(x,y)], pch=0)
points(l[type=='s',.(x,y)], pch='N', col='red', cex=2)

# multi-panel version
plotter_by_risk = function(base_filename, risk_threshold, pnet, as.png=T) {
  snet = pnet[r > risk_threshold]
  plotter(paste0(base_filename, risk_threshold), snet, 'black', as.png)
}

png('interaction_map_5-panels.png', width=10000, height=2000, res=500)
par(mfrow=c(1,5))
for (ppb in seq(0,1,0.25)) {
  plotter_by_risk('', ppb, plot_net, as.png=F)
}
dev.off()

png('interaction_map_3-panels.png', width=6000, height=2000, res=400)
par(mfrow=c(1,3))
for (ppb in c(0, 0.5, 0.75)) {
  plotter_by_risk('', ppb, plot_net, as.png=F)
}
dev.off()

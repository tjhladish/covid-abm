require('rjson')
require("RSQLite")

vers=NA
arg=commandArgs(trailingOnly=T)
if (length(arg) == 0) {
  cat("Argument not found. Using default version number and number of particles.\n")
  numparticle=500
} else if (length(arg) != 2) {
  stop("Must have 2 arguments: version number and number of particles.\n")
} else {
  vers=arg[1]
  numparticle=as.numeric(arg[2])
}

# parse json
abc_config=fromJSON(file="abc_covid.json");
if (is.na(vers)) vers=sub(".*_v(.*).sqlite", "\\1", abc_config$database_filename)
nmet=length(abc_config$metrics);
npar=length(abc_config$parameters);
#factors=c(rep(1,nmet), rep(2,npar)); factors %o% factors
abc_metrics_values=rep(0,length(abc_config$metrics)); for(i in 1:length(abc_config$metrics)) { abc_metrics_values[i]=abc_config$metrics[[i]]$value; }
abc_metrics_names=rep(0,length(abc_config$metrics)); for(i in 1:length(abc_config$metrics)) { abc_metrics_names[i]=abc_config$metrics[[i]]$name; }
# we want short names if they're defined
for(i in 1:length(abc_config$metrics)) { if('short_name' %in% names(abc_config$metrics[[i]])) abc_metrics_names[i]=abc_config$metrics[[i]]$short_name; }

# read in db
drv = dbDriver("SQLite")
db = dbConnect(drv, abc_config$database_filename)
abc = dbGetQuery(db, 'select J.*, P.*, M.* from job J, par P, met M where J.serial = P.serial and J.serial = M.serial and smcSet = (select max(smcSet) from job) and posterior > -1')
#abc = dbGetQuery(db, 'select J.*, P.*, M.* from job J, upar P, met M where J.serial = P.serial and J.serial = M.serial and smcSet = 6 and posterior > -1')
extra_serials = which(names(abc)== 'serial')[-1]
abc = abc[,-c(extra_serials)]
abc = subset(abc, select=-c(startTime, duration, attempts, seed))

col_offset = 5
par_cols = 1:npar + col_offset
met_cols = 1:nmet + col_offset + npar

#abc$above = abc$deaths>9+4*abc$T
#abc$above = 1 + (abc$rcases_all > 5000)
abc$color = colorRampPalette(c("green", "grey", "black"))( numparticle )[abc$posterior + 1]

proto = abc[1,]
proto[1,] <- NA
proto[1,abc_metrics_names] <- abc_metrics_values
dm = rbind(abc, proto)
dm$sim = factor(ifelse(is.na(dm$serial), F, T))
dm[dm$sim==F, par_cols] = apply(dm[dm$sim==T,par_cols], 2, median)

#colors <- c(sim = '#00000012', par = 'purple', met = 'orange')
#chars <- c(sim = 20, par = '|', met = '—')
############## end



#names(dm)[4:12] = c('EF', 'Mos move', 'Daily intros', 'Num mos', 'Mean', 'Median', 'Stdev', 'Max', 'Skew')
#names(dm)[4:13] = c('EF', 'Mos move', 'Daily intros', 'Num mos', 'Beta', 'Mean', 'Median', 'Stdev', 'Max', 'Skew')
#names(dm)[13] = "Autocorr"
#names(dm)[14] = "Autocorr"


source("pairs.panels.R")
#pdf("pairs-a.pdf", width=16, height=16)
outfile=paste0("pairs-par_v", vers, ".png")
png(outfile, width=5000, height=3000, res=240)
pairs.panels(dm[,c(par_cols, met_cols)], dm[,dim(dm)[2]], npar, nmet, points.col=dm$color,# '#00000020',
             box.col='black', box.lwd=0.5, gap=0.5, cor=F, points.cex=0.5)
dev.off()

#miniplot_subset =c(par_cols[c(5,6,3)], met_cols[c(4,11)])
# names(dm)[miniplot_subset] = c('Mosq/loc', 'Bite & Transmit', "Post-pri. severity", 'Median','Severe prev.')
# png("pairs-b.png", width=1800, height=1340, res=180)
# pairs.panels(dm[,miniplot_subset], dm[,dim(dm)[2]], npar=3, nmet=2, points.col='#00000012',
#              box.col='black', box.lwd=0.5, gap=0.5, cor=F, line_wt=2)
# dev.off()

## UPAR
abc = dbGetQuery(db, 'select J.*, P.*, M.* from job J, upar P, met M where J.serial = P.serial and J.serial = M.serial and smcSet = 6 and posterior > -1')
extra_serials = which(names(abc)== 'serial')[-1]
abc = abc[,-c(extra_serials)]
abc = subset(abc, select=-c(startTime, duration, attempts, seed))

col_offset = 5
par_cols = 1:npar + col_offset
met_cols = 1:nmet + col_offset + npar

#abc$above = abc$deaths>9+4*abc$T
#abc$above = 1 + (abc$rcases_all > 5000)
abc$color = colorRampPalette(c("green", "grey", "black"))( numparticle )[abc$posterior + 1]

proto = abc[1,]
proto[1,] <- NA
proto[1,abc_metrics_names] <- abc_metrics_values
dm = rbind(abc, proto)
dm$sim = factor(ifelse(is.na(dm$serial), F, T))
dm[dm$sim==F, par_cols] = apply(dm[dm$sim==T,par_cols], 2, median)

#colors <- c(sim = '#00000012', par = 'purple', met = 'orange')
#chars <- c(sim = 20, par = '|', met = '—')
############## end



#names(dm)[4:12] = c('EF', 'Mos move', 'Daily intros', 'Num mos', 'Mean', 'Median', 'Stdev', 'Max', 'Skew')
#names(dm)[4:13] = c('EF', 'Mos move', 'Daily intros', 'Num mos', 'Beta', 'Mean', 'Median', 'Stdev', 'Max', 'Skew')
#names(dm)[13] = "Autocorr"
#names(dm)[14] = "Autocorr"


source("pairs.panels.R")
#pdf("pairs-a.pdf", width=16, height=16)
outfile=paste0("pairs-upar_v", vers, ".png")
png(outfile, width=5000, height=3000, res=240)
pairs.panels(dm[,c(par_cols, met_cols)], dm[,dim(dm)[2]], npar, nmet, points.col=dm$color,# '#00000020',
             box.col='black', box.lwd=0.5, gap=0.5, cor=F, points.cex=0.5)
dev.off()

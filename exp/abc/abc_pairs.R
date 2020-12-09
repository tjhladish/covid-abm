require('rjson')
require("RSQLite")

json_file = 'abc_config.json'
vers=NA
arg=commandArgs(trailingOnly=T)
if (length(arg) == 0) {
    cat("Arguments not provided. Using 'abc_config.json' and default version number if available.\n")
} else if (length(arg) != 2) {
    stop("Expecting have 0 or 2 arguments: JSON filename and version number.\n")
} else {
    json_file = arg[1]
    vers=arg[2]
}

# parse json
abc_config=fromJSON(file=json_file);
if (is.na(vers)) vers=sub(".*_v(.*).sqlite", "\\1", abc_config$database_filename)

nmet=length(abc_config$metrics);
npar=length(abc_config$parameters);
col_offset = 5
par_cols = 1:npar + col_offset
met_cols = 1:nmet + col_offset + npar

#factors=c(rep(1,nmet), rep(2,npar)); factors %o% factors
abc_metrics_values=rep(0,length(abc_config$metrics)); for(i in 1:length(abc_config$metrics)) { abc_metrics_values[i]=abc_config$metrics[[i]]$value; }
abc_metrics_names=rep(0,length(abc_config$metrics)); for(i in 1:length(abc_config$metrics)) { abc_metrics_names[i]=abc_config$metrics[[i]]$name; }
# we want short names if they're defined
for(i in 1:length(abc_config$metrics)) { if('short_name' %in% names(abc_config$metrics[[i]])) abc_metrics_names[i]=abc_config$metrics[[i]]$short_name; }

# read in db
drv = dbDriver("SQLite")
db = dbConnect(drv, abc_config$database_filename)
abc_par = dbGetQuery(db, 'select J.*, P.*, M.* from job J, par P, met M where J.serial = P.serial and J.serial = M.serial and smcSet = (select max(smcSet) from job where posterior > -1) and posterior > -1')
abc_upar = dbGetQuery(db, 'select J.*, P.*, M.* from job J, upar P, met M where J.serial = P.serial and J.serial = M.serial and smcSet = (select max(smcSet) from job where posterior > -1) and posterior > -1')
dbDisconnect(db)

numparticle = nrow(abc_par)
extra_serials = which(names(abc_par)== 'serial')[-1]

source("pairs.panels.R")

plotter = function(abc, fileroot) {
    abc = abc[,-c(extra_serials)]
    abc = subset(abc, select=-c(startTime, duration, attempts, seed))

    abc$color = colorRampPalette(c("black", "grey", "lightgrey"))( numparticle )[abc$posterior + 1]

    proto = abc[1,]
    proto[1,] <- NA
    proto[1,abc_metrics_names] <- abc_metrics_values
    dm = rbind(abc, proto)
    dm$sim = factor(ifelse(is.na(dm$serial), F, T))
    dm[dm$sim==F, par_cols] = apply(dm[dm$sim==T,par_cols], 2, median)

    #outfile=paste0("pairs-par_v", vers, ".png")
    outfile=paste0(fileroot, vers, ".png")
    png(outfile, width=5000, height=3000, res=240)
    pairs.panels(dm[,c(par_cols, met_cols)], dm[,dim(dm)[2]], npar, nmet, points.col=dm$color,# '#00000020',
                 box.col='black', box.lwd=0.5, gap=0.5, cor=F, points.cex=0.5)
    dev.off()
}

plotter(abc_par, 'pairs-par_v')
plotter(abc_upar, 'pairs-upar_v')

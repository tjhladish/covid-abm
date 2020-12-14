require('rjson')
require("RSQLite")
require(beanplot)

json_file = 'abc_covid.json'
vers=NA
arg=commandArgs(trailingOnly=T)
if (length(arg) == 0) {
    cat("Arguments not provided. Using 'abc_covid.json' and default version number if available.\n")
} else if (length(arg) == 1) {
	json_file = arg[1]
	cat("Only first argument is provided. Inferring version number from the .json file.\n")
} else if (length(arg) == 2) {
    json_file = arg[1]
    vers=arg[2]
} else {
	stop("Expecting 0, 1 or 2 arguments: JSON filename and version number.\n")
}

# Functions
preprocess = function (abc) {
	extra_serials = which(names(abc)== 'serial')[-1]
	abc = abc[,-c(extra_serials)]
	abc = subset(abc, select=-c(startTime, duration, attempts, seed))
	abc$post_bool = abc$posterior >= 0

	col_offset = 5
	par_cols = 1:npar + col_offset
	met_cols = 1:nmet + col_offset + npar
	obs_met  = sapply(abc_config$metrics, '[[', 'value')#c(16.7, 322, 349, 5.26, 201, 0.160, -0.106, 1.23, 6.02, 9.75, 21.7, 0.98, 1.00, 1.00, 0.78, 1.37)

	incomplete_sets = unique(abc$smcSet[abc$status!='D']) # 6
	all_sets = unique(abc$smcSet)
	complete_sets = setdiff(all_sets, incomplete_sets)
	last_complete_set = max(complete_sets)
	
	return(list(abc = abc, par_cols = par_cols, 
				met_cols = met_cols, obs_met = obs_met,
				last_complete_set = last_complete_set))
}

plotter = function (abc, type = "par", vers) {
	outfile = paste0('marginal_', type, 's_v', vers, '.pdf')
	pdf(outfile, width=8.5, height=11)
	par(mfrow=c(npar,1))
	par(mar=c(2.1, 4.1, 1.1, 0.5))
	for (col in par_cols) {
		colname = names(abc)[col];
		cat(paste(colname, '\n'))
		beanplot( abc[,colname] ~ abc$post_bool*abc$smcSet, what=c(1,1,1,0), col=list(c(grey(0.3), 1,1, grey(0)), c(grey(0.9), 1,1, grey(0.7))), main='', ylab=colname, side='both')
	}
	dev.off()
}

# parse json
abc_config=fromJSON(file=json_file);
if (is.na(vers)) vers=sub(".*_v(.*).sqlite", "\\1", abc_config$database_filename)
nmet=length(abc_config$metrics);
npar=length(abc_config$parameters);

drv = dbDriver("SQLite")
db = dbConnect(drv, abc_config$database_filename)
tables = dbListTables(db)
abc = dbGetQuery(db, 'select J.*, P.*, M.* from job J, par P, met M where J.serial = P.serial and J.serial = M.serial')
tmp = preprocess(abc)
abc = tmp$abc
par_cols = tmp$par_cols
met_cols = tmp$met_cols
obs_met = tmp$obs_met
last_complete_set = tmp$last_complete_set

if ("upar" %in% tables) {
	abc_upar = dbGetQuery(db, 'select J.*, P.*, M.* from job J, upar P, met M where J.serial = P.serial and J.serial = M.serial')
	abc_upar = preprocess(abc_upar)$abc
}
dbDisconnect(db)

# plot
plotter(abc, vers = vers)
if (exists("abc_upar")) plotter(abc_upar, type = "upar", vers = vers)


#ylims = list(c(0,800), c(0,5), c(0,200), c(0,350), c(0,600), c(10,30000), c(10e-1, 10e4), c(-2, 6), c(0, 1), c(0, 1), c(0, 1), c(-6, 4), c(-0.1, 0.1))

# 100 for these values indicates an inability to do a logistic fit.  This is rare, but could be avoided if we made the regression approach more robust.
#abc$beta0[abc$beta0 == 100] = NA
#abc$beta1[abc$beta1 == 100] = NA

# Plot metrics - tends to be trickier, because distributions can be weird (e.g, highly skewed or long-tailed)
outfile = paste0('marginal_mets_v', vers, '.pdf')
pdf(outfile, width=8.5, height=50)
par(mfrow=c(nmet,1))
par(mar=c(2.1, 4.1, 1.1, 0.5))
alpha = 0.025 # plot middle 95% of distributions
for (col in met_cols) {
    met_idx = col - max(par_cols)
    colname = names(abc)[col]

    # calculate reasonable plot limits
    obs_val = obs_met[met_idx]
    val_lims = quantile(abc[,colname], na.rm=T, probs=c(alpha, 1-alpha))
    val_min = min(val_lims[1], obs_val)
    val_max = max(val_lims[2], obs_val)

    # filter out NAs
    complete = complete.cases(abc[,colname]) & abc[,colname] > val_lims[1] & abc[,colname] < val_lims[2]
    cat(paste0(colname, ' ', val_lims[1], ', ', val_lims[2], '\n'))

    # plot the stuff
    beanplot( abc[complete, colname] ~ abc$post_bool[complete] * abc$smcSet[complete], what=c(0,1,1,0),
              col=list(c(grey(0.3), 1,1, grey(0)), c(grey(0.9), 1, 1, grey(0.7))),
              #main='', ylab=colname, side='both', ylim=unlist(ylims[met_idx]) )
              main='', ylab=colname, side='both', ylim=c(val_min, val_max), na.rm=T, log='' )
    abline(h=mean(abc[abc$smcSet==last_complete_set, colname], na.rm=T), lty=3)
    abline(h=obs_val, col=2, lwd=1)
}

dev.off()

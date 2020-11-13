library(data.table)

# Script for processing SVG line data grabbed from https://covid19.gleamproject.org/mobility

# gleam_cmi-raw.csv    Cuebiq mobility index
# gleam_cuci-raw.csv   Cuebiq unique contact index
# gleam_ccdi-raw.csv   Cuebiq contact duration index

DATE_FIRST = as.Date('2020-02-01') # first day in timeseries
DATE_LAST  = as.Date('2020-10-25') # last day in timeseries

process_gleam = function(filename, metric_name, metric_min, metric_max, date_min, date_max) {
    # read file, convert to 2-col data table
    d = as.numeric(readLines(filename))
    m = matrix(d, ncol=2, byrow=T)
    dt = data.table(m)
    names(dt) = c('x', 'y')

    # coordinate space from the SVG object is a linear transformation of the actual data
    # therefore, we need to re-normalize:
    y_min = min(dt$y)
    y_max = max(dt$y)
    dt[,metric_name] = (1 - (dt$y - y_min)/(y_max-y_min))*(metric_max - metric_min) + metric_min

    x_min = min(dt$x)
    x_max = max(dt$x)
    dt$date = (dt$x - x_min)/(x_max-x_min)*(date_max - date_min) + date_min

    # dt$date ends up being a bit wonky, because of daylight savings time and the fact that POSIX knows about it
    stopifnot(nrow(dt) %% 3 == 1) # We're expecting three values for each day, except the last one, which we expect to have one
    dt$date2 = c(rep(seq(DATE_FIRST, DATE_LAST - 1, by="day"), each = 3), DATE_LAST)

    # average over the values reported (typically 3 of them) for each day
    model_daily = dt[, .(metric_name = mean(get(metric_name))), by=.(date=date2)]
    max_rel_metric = max(model_daily[which(model_daily$date >= as.Date('2020-03-01')),2])

    adj_name = paste0('adj_', metric_name)
    setnames(model_daily, "metric_name", metric_name)  # "metric_name" --> e.g. "cuci"

    # rescale metric to have a max of 1, as that is what the simulator expects
    model_daily[,adj_name] = model_daily[,metric_name, with=F]/max_rel_metric # this "with" stuff is some serious bullshit
    print(model_daily)

    # remove column for unscaled version
    model_daily[,(metric_name):=NULL]

    # set early values to their mean, so that we aren't dealing with extreme sensitivity to starting date
    model_daily[which(model_daily$date < as.Date('2020-03-10')), adj_name] = mean(as.numeric(unlist(model_daily[date < as.Date('2020-03-10'), adj_name, with=F])))

    return(model_daily);
}

date_min = as.POSIXct("2020-02-01")  # first date reported on GLEAM website
date_max = as.POSIXct("2020-10-25")  # last date reported on GLEAM website

# min and max values read from tooltip on https://covid19.gleamproject.org/mobility
cmi_dt  = process_gleam('gleam_cmi-raw.csv',  'cmi',  metric_min = 0.255, metric_max = 1.315, date_min, date_max) # mobility index
cuci_dt = process_gleam('gleam_cuci-raw.csv', 'cuci', metric_min = 0.387, metric_max = 3.314, date_min, date_max) # unique contact index
ccdi_dt = process_gleam('gleam_ccdi-raw.csv', 'ccdi', metric_min = 0.314, metric_max = 3.431, date_min, date_max) # contact duration index

#cci_min = 0.1765 # values for number of contacts metric
#cci_max = 1.7658
#cci_min = 0.1654
#cci_max = 1.3195

model = data.table(date = seq(as.Date("2020-01-01"), DATE_LAST, by="day"), compound = as.numeric(NA))
model = merge(model, cmi_dt,  by='date', all.x=T)
model = merge(model, cuci_dt, by='date', all.x=T)
model = merge(model, ccdi_dt, by='date', all.x=T)
model$compound = (model$adj_cmi * model$adj_cuci * model$adj_ccdi)**(1/3)
#model$compound = model$compound/max(model$compound, na.rm=T)
nc = ncol(model)
model[is.na(model$compound), 2:nc] = model[model$date == as.Date('2020-02-01'), 2:nc]

write.csv(model, "gleam_metrics.csv", quote=F, row.names=F)

png("gleam_metrics.png", width = 2400, height = 1600, res=200)

plot (model$date, model$compound, type='l', lwd=3, ylim=c(0,1), xlab='Date', ylab='Index value')
lines(model$date, model$adj_ccdi, col='blue')
lines(model$date, model$adj_cuci, col='red')
lines(model$date, model$adj_cmi, col='green')
legend("bottomright", bty='n', legend=c('Cuebiq unique contact index', 'Cuebiq contact duration index', 'Cuebiq mobility index', 'Cube-root of product'), col=c('red', 'blue', 'green', 'black'), lwd=3)
dev.off()

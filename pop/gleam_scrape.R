library(data.table)

# Script for processing SVG line data grabbed from https://covid19.gleamproject.org/mobility

# gleam_cmi-raw.csv    Cuebiq mobility index
# gleam_cuci-raw.csv   Cuebiq unique contact index
# gleam_ccdi-raw.csv   Cuebiq contact duration index

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

    adj_name = paste0('adj_', metric_name)
    model_daily = dt[, .(metric_name = mean(get(metric_name))), by=.(date=date2)]
    setnames(model_daily, "metric_name", metric_name)  # "metric_name" --> e.g. "cuci"

    # set early values to their mean, so that we aren't dealing with extreme sensitivity to starting date
    model_daily[which(model_daily$date < as.Date('2020-03-10')), metric_name] = mean(as.numeric(unlist(model_daily[date < as.Date('2020-03-10'), metric_name, with=F])))

    # average over the values reported (typically 3 of them) for each day
    max_rel_metric = max(model_daily[which(model_daily$date >= as.Date('2020-03-01')),2])

    # rescale metric to have a max of 1, as that is what the simulator expects
    model_daily[,adj_name] = 1 - model_daily[,metric_name, with=F]/max_rel_metric # this "with" stuff is some serious bullshit
    print(model_daily)

    # remove column for unscaled version
    model_daily[,(metric_name):=NULL]

    return(model_daily);
}

DATE_FIRST = as.Date('2020-02-01') # first day in timeseries
DATE_LAST  = as.Date('2021-03-04') # last day in timeseries

date_min = as.POSIXct("2020-02-01")  # first date reported on GLEAM website
date_max = as.POSIXct("2021-03-04")  # last date reported on GLEAM website

# min and max values read from tooltip on https://covid19.gleamproject.org/mobility
cuci_dt = process_gleam('gleam_scrape-cuci-raw-210501.txt', 'cuci', metric_min = 0.387, metric_max = 3.314, date_min, date_max) # unique contact index

model = data.table(date = seq(as.Date("2020-01-01"), DATE_LAST, by="day"))
model = merge(model, cuci_dt, by='date', all.x=T)
model$adj_cuci[is.na(model$adj_cuci)] = model$adj_cuci[model$date == as.Date('2020-02-01')]
nc = ncol(model)

write.csv(model, "gleam_metrics.csv", quote=F, row.names=F)

png("gleam_metrics.png", width = 2400, height = 1600, res=200)

plot (model$date, model$adj_cuci, type='l', lwd=3, ylim=c(0,1), xlab='Date', ylab='Index value')
legend("bottomright", bty='n', legend=c('Cuebiq unique contact index'), lwd=3)
dev.off()

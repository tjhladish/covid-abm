library(data.table)

# gleam_cmi-raw.csv    Cuebiq mobility index
# gleam_cuci-raw.csv   Cuebiq unique contact index
# gleam_ccdi-raw.csv   Cuebiq contact duration index

DATE_FIRST = as.Date('2020-02-01') # first day in timeseries
DATE_LAST  = as.Date('2020-10-25') # last day in timeseries

process_gleam = function(filename, metric_name, metric_min, metric_max, date_min, date_max) {
    d = as.numeric(readLines(filename))
    m = matrix(d, ncol=2, byrow=T)
    dt = data.table(m)
    names(dt) = c('x', 'y')
    #cci_min = 0.1654
    #cci_max = 1.3195
    y_min = min(dt$y)
    y_max = max(dt$y)
    dt[,metric_name] = (1 - (dt$y - y_min)/(y_max-y_min))*(metric_max - metric_min) + metric_min

    x_min = min(dt$x)
    x_max = max(dt$x)
    dt$date = (dt$x - x_min)/(x_max-x_min)*(date_max - date_min) + date_min

    stopifnot(nrow(dt) %% 3 == 1) # We're expecting three values for each day, except the last one, which we expect to have one
    dt$date2 = c(rep(seq(DATE_FIRST, DATE_LAST - 1, by="day"), each = 3), DATE_LAST)
    model_daily = dt[, .(metric_name = mean(get(metric_name))), by=.(date=date2)]
    max_rel_metric = max(model_daily[which(model_daily$date >= as.Date('2020-03-01')),2])

    adj_name = paste0('adj_', metric_name)
    setnames(model_daily, "metric_name", metric_name)
    model_daily[,adj_name] = model_daily[,metric_name, with=F]/max_rel_metric # this "with" stuff is some serious bullshit
    print(model_daily)
    model_daily[,(metric_name):=NULL]

    return(model_daily);
}

date_min = as.POSIXct("2020-02-01")
date_max = as.POSIXct("2020-10-25")

cmi_dt  = process_gleam('gleam_cmi-raw.csv',  'cmi',  metric_min = 0.255, metric_max = 1.315, date_min, date_max)
cuci_dt = process_gleam('gleam_cuci-raw.csv', 'cuci', metric_min = 0.387, metric_max = 3.314, date_min, date_max)
ccdi_dt = process_gleam('gleam_ccdi-raw.csv', 'ccdi', metric_min = 0.314, metric_max = 3.431, date_min, date_max)

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
model$compound[is.na(model$compound)] = model$compound[model$date == as.Date('2020-02-01')]

write.csv(model, "gleam_metrics.csv", quote=F, row.names=F)

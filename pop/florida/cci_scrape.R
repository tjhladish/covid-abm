library(data.table)

d = as.numeric(readLines('gleam_cci-raw.csv'))
m = matrix(d, ncol=2, byrow=T)
cci_dt = data.table(m)
names(cci_dt) = c('x', 'y')
cci_min = 0.1765
cci_max = 1.7658
y_min = min(cci_dt$y)
y_max = max(cci_dt$y)
cci_dt$cci = (1 - (cci_dt$y - y_min)/(y_max-y_min))*(cci_max - cci_min) + cci_min

d2 = as.numeric(readLines('gleam_cmi-raw.csv'))
m2 = matrix(d2, ncol=2, byrow=T)
cmi_dt = data.table(m2)
names(cmi_dt) = c('x', 'y')
cmi_min = 0.2488 
cmi_max = 1.2985
y_min = min(cmi_dt$y)
y_max = max(cmi_dt$y)
cmi_dt$cmi = (1 - (cmi_dt$y - y_min)/(y_max-y_min))*(cmi_max - cmi_min) + cmi_min

date_min = as.POSIXct("2020-02-01")
date_max = as.POSIXct("2020-10-25")

x_min = min(cci_dt$x)
x_max = max(cci_dt$x)
cci_dt$date = (cci_dt$x - x_min)/(x_max-x_min)*(date_max - date_min) + date_min
cci_dt$date2 = c(rep(seq(as.Date("2020-02-01"), as.Date("2020-10-24"), by="day"), each = 3), as.Date("2020-10-25"))
model_cci = cci_dt[, .(cci = mean(cci)), by=.(date=date2)]
max_rel_cci = max(model_cci[which(model_cci$date >= as.Date('2020-03-01')),2])
model_cci$adjusted_cci = model_cci$cci/max_rel_cci
#model_cci$adjusted_cci[which(model_cci$date < as.Date('2020-03-10'))] = 1.0
model_cci$adjusted_cci[which(model_cci$date < as.Date('2020-03-10'))] = as.numeric(model_cci[date == as.Date('2020-03-10'),'adjusted_cci'])

x_min = min(cmi_dt$x)
x_max = max(cmi_dt$x)
cmi_dt$date = (cmi_dt$x - x_min)/(x_max-x_min)*(date_max - date_min) + date_min
cmi_dt$date2 = c(rep(seq(as.Date("2020-02-01"), as.Date("2020-10-24"), by="day"), each = 3), as.Date("2020-10-25"))
model_cmi = cmi_dt[, .(cmi = mean(cmi)), by=.(date=date2)]
max_rel_cmi = max(model_cmi[which(model_cmi$date >= as.Date('2020-03-01')),2])
model_cmi$adjusted_cmi = model_cmi$cmi/max_rel_cmi
model_cmi$adjusted_cmi[which(model_cmi$date < as.Date('2020-03-10'))] = as.numeric(model_cmi[date == as.Date('2020-03-10'),'adjusted_cmi'])

model = data.table(date = seq(as.Date("2020-01-01"), as.Date("2020-10-25"), by="day"), rs_cci_cmi = as.numeric(NA))
model = merge(model, model_cci, by='date', all.x=T)
model = merge(model, model_cmi, by='date', all.x=T)
model$rs_cci_cmi = sqrt(model$adjusted_cci * model$adjusted_cmi)
model$rs_cci_cmi = model$rs_cci_cmi/max(model$rs_cci_cmi, na.rm=T)
model$rs_cci_cmi[is.na(model$rs_cci_cmi)] = 1.0

#write.csv(model[,c(1,3,2)], "gleam_cci_cmi-adjusted.csv", quote=F, row.names=F)

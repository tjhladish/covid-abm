library('data.table')
library('lubridate')

d = fread('Daily_Sero_with_Type_data.csv')
d2 = data.table()
d2$date = lubridate::mdy(d$'Day of Incidence Date')
d2$lower = d$'Avg. Cumulative Prevalence Lower CI'/100.0
d2$est = d$'Avg. Cumulative Prevalence Rate (%)'/100.0
d2$upper = d$'Avg. Cumulative Prevalence Upper CI'/100.0

write.csv(data.frame(d2), "CDC_seroprev_long.csv", quote=F, row.names=F)

if (interactive()) { setwd("~/documents/work/covid-abm/exp/active-vac/") }

.args <- if (interactive()) c(
  "rcasedeath-florida.csv"
) else commandArgs(trailingOnly = TRUE)

ed = read.csv(.args[1])
ed$date = as.Date(ed$Date)

ed$crcase  = cumsum(ed$rcase)
ed$crdeath = cumsum(ed$rdeath)

ed$rcase_ma7 = filter(c(0,0,0,ed$rcase), rep(1/7, 7), sides = 2)[-(1:3)]
ed$crcase_ma7 = cumsum(ed$rcase_ma7)

out = data.frame(Date = ed$Date, rcase_ma7 = ed$rcase_ma7, rdeath = ed$rdeath, death_incd = ed$death_incd, excess = ed$excess)
out = out[1:(nrow(out)-3),]
write.csv(out, "rcasedeath-florida_ma7cases.csv", quote = F, row.names = F)

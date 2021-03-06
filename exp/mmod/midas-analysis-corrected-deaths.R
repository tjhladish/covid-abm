importData = read.csv("mmod_v2.3.csv", header=T, sep='\t', stringsAsFactors = F)

#0 = "closed"
#1 = "2weeks"
#2 = "5percent"
#3 = "open"

anyOutbreak <- function(x){
  return(as.numeric(any(x >= 10)))
}

days_closed <- aggregate(importData$closed, by = list(scen = importData$scen, rep = importData$rep), sum)
cumu_infections <- aggregate(importData$infinc, by = list(scen = importData$scen, rep = importData$rep), sum)
cumu_deaths <- aggregate(importData$rdeath, by = list(scen = importData$scen, rep = importData$rep), sum)
peak_hosp <- aggregate(importData$hosprev, by = list(scen = importData$scen, rep = importData$rep), max)
prob_outbreak <- aggregate(importData$rcases, by = list(scen = importData$scen, rep = importData$rep), anyOutbreak)

daysQuantile <- aggregate(days_closed$x, by = list(days_closed$scen), FUN = 'quantile', probs = seq(0, 1, 0.01))
infectionsQuantile <- aggregate(cumu_infections$x, by = list(cumu_infections$scen), FUN = 'quantile', probs = seq(0, 1, 0.01))
deathsQuantile <- aggregate(cumu_deaths$x, by = list(cumu_deaths$scen), FUN = 'quantile', probs = seq(0, 1, 0.01))
hospQuantile <- aggregate(peak_hosp$x, by = list(peak_hosp$scen), FUN = 'quantile', probs = seq(0, 1, 0.01))
# outbreakQuantile <- aggregate(prob_outbreak$x, by = list(prob_outbreak$scen), FUN = 'quantile', probs = seq(0, 1, 0.01))
outbreakProbs <- aggregate(prob_outbreak$x, by = list(prob_outbreak$scen), FUN = 'mean')

#row 1 = "closed", row 2 = "onepercent", row 3 = open", row 4 = "twoweeks"
closedQuantiles <- data.frame(quantile = 0:100, intervention = "closed", cumu_infections = infectionsQuantile$x[infectionsQuantile$Group.1 == 0], 
                              cumu_deaths = deathsQuantile$x[deathsQuantile$Group.1 == 0], peak_hosp = hospQuantile$x[hospQuantile$Group.1 == 0], 
                              prob_outbreak = NA, days_closed = daysQuantile$x[daysQuantile$Group.1 == 0])
closedQuantiles$prob_outbreak[closedQuantiles$quantile == 50] <- outbreakProbs$x[outbreakProbs$Group.1 == 0]

twoWeeksQuantiles <- data.frame(quantile = 0:100, intervention = "2weeks", cumu_infections = infectionsQuantile$x[infectionsQuantile$Group.1 == 1], 
                                cumu_deaths = deathsQuantile$x[deathsQuantile$Group.1 == 1], peak_hosp = hospQuantile$x[hospQuantile$Group.1 == 1], 
                                prob_outbreak = NA, days_closed = daysQuantile$x[daysQuantile$Group.1 == 1])
twoWeeksQuantiles$prob_outbreak[twoWeeksQuantiles$quantile == 50] <- outbreakProbs$x[outbreakProbs$Group.1 == 1]

onePercentQuantiles <- data.frame(quantile = 0:100, intervention = "5percent", cumu_infections = infectionsQuantile$x[infectionsQuantile$Group.1 == 2], 
                                  cumu_deaths = deathsQuantile$x[deathsQuantile$Group.1 == 2], peak_hosp = hospQuantile$x[hospQuantile$Group.1 == 2], 
                                  prob_outbreak = NA, days_closed = daysQuantile$x[daysQuantile$Group.1 == 2])
onePercentQuantiles$prob_outbreak[onePercentQuantiles$quantile == 50] <- outbreakProbs$x[outbreakProbs$Group.1 == 2]

openQuantiles <- data.frame(quantile = 0:100, intervention = "open", cumu_infections = infectionsQuantile$x[infectionsQuantile$Group.1 == 3], 
                            cumu_deaths = deathsQuantile$x[deathsQuantile$Group.1 == 3], peak_hosp = hospQuantile$x[hospQuantile$Group.1 == 3], 
                            prob_outbreak = NA, days_closed = daysQuantile$x[daysQuantile$Group.1 == 3])
openQuantiles$prob_outbreak[openQuantiles$quantile == 50] <- outbreakProbs$x[outbreakProbs$Group.1 == 3]

analysisOutput <- data.frame(team = "UF", round = 2, rbind(closedQuantiles, twoWeeksQuantiles, onePercentQuantiles, openQuantiles))

write.csv(analysisOutput, file = "MIDAS-v2.3.csv", row.names = FALSE)

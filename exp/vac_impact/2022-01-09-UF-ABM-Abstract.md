# Summary of results
Using an agent-based model with demographics, vaccination coverage, and pandemic history for the state of Florida, we projected the cases, hospitalizations, and deaths caused by the SARS-CoV-2 Omicron variant.
Our model ranks the proposed scenarios, from most to least optimistic, in the following order: SCENARIO RANK
BRIEF DESCRIPTION OF OUTCOMES FOR DIFFERENT SCENARIOS

# Explanation of observed dynamics given model assumptions
One important assumption in our model is that fully vaccinated individuals are less infectious compared to partially vaccinated and unvaccinated individuals.
Given that approximately 60% of the population is fully vaccinated prior to the Omicron wave, those with breakthrough Omicron infections will experience less onward transmission.
Thus, in THIS STILL MAY HOLD TRUE FOR CERTAIN SCENARIOS, given a higher immune escape assumption but lower overall transmission advantage, even though many breakthrough Omicron infections may occur, these will result in fewer secondary infections than in infections in unvaccinated people, partially explaining the more optimistic projections for these scenarios.

# Model assumptions
## Initial distribution of susceptibility
Initial distribution of susceptibility is constant across all replicates up until December 1, 2021, when we assume Omicron first started to be introduced to the state.
We explicitly model all previous waves of the pandemic, and a vaccine-roll out based on what has occurred in Florida (i.e., health care workers first, followed age-eligibility changes as they occurred in the state).

## Transmillibily

## Generation time
Our model assumes that the wildtype virus and non-Omicron mutants have average latent and incubation periods of XX and XX days (respectively).
This translates to an average modeled generation time of XX days.
For the Omicron VOC, we assume the average latent and incubation periods are reduced to XX and XX days, resulting in an average generation time of XX days.

## Waning immunity assumptions
In past work, we determined that assuming leaky immunity and changes in VOC transmissibility and immune escape ability were better at predicting apparent loss of vaccine protection than immunity that intrinsically wanes.
Our model therefore does not assume intrinsic waning of immunity from prior infection or vaccination.
Our assumptions result in estimated vaccine direct effectiveness against infection that is ~0.8 early in 2021, down to ~0.2 against delta and ~0.1 against Omicron.

## Seasonality
Seasonality is modeled using a sine model with a 6 month period, with transmission probability peaking in July and January.

## Baseline conditions for delta VOC prevalence
We assume that when Omicron is first introduced in our model (December 1, 2021), infections are very low, but all are caused by delta.
This is consistent with Florida data.

## Changing severity with VOCs
We interpret severity as the probability that a symptomatic infection will result in a severe (i.e., hospitalizable) infection.
We assume alpha severity is similar to wildtype infections, whereas delta severity is 2.5x higher.
Omicron severity was defined relative to delta, according to the specified scenarios.

## Implementing NPIs
Closure of non-essential businesses and schools occur in our model early in the pandemic, but not during the Omicron wave.
Individuals in the model have time-varying personal-protective behaviors, namely avoiding social contacts and high-risk businesses like bars and restaurants.

# Other updates in model assumptions from previous rounds
## Changes in reporting outcomes due to Omicron
Infections become detected with a probability dependent on the severity of symptoms and which changes over time.
The chance of detection at each state (asymptomatic, mild, severe, critical, or dead) is represented as a probability of detection at that state if not detected at previous states.
Below these probabilities are listed for given time points in the simulation (the model extrapolates between these points). 

| Date       | Asymptomatic | Mild  | Severe | Critical | Dead     |
| :---:      | :---:        | :---: | :---:  | :---:    | :---:    |
| 2020-02-05 | 0.0          | 0.05  | 0.7    | 0.1      | 0.220273 |
| 2020-06-01 | 0.02         | 0.5   | 0.1    | 0.1      | 0.496095 |
| 2020-10-01 | 0.1          | 0.9   | 0.9    | 0.0      | 1.0      |
| 2021-11-01 | 0.02         | 0.9   | 0.9    | 0.0      | 1.0      |

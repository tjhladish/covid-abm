# covid-abm
---
This stochastic, discrete-time, agent-based model of dengue transmission explicitly represents a synthetic population of humans with which transmission of SARS-CoV-2 (including variants of concern) can be simulated.

This model is adapted from a previous [agent-based model for dengue transmission](https://github.com/tjhladish/dengue).

## Daily Transmission Cycle
The model revolves around a daily cycle of human movement and transmission opportunities.
Susceptible humans can become infected when they co-localize with currently infectious individuals.
Infections follow an SEIRD model progressing between **S**usceptible, **E**xposed, **I**nfectious, and **R**ecovered, and **D**ead states. 
Moreover, infections can develop varying levels of severity: asymptomatic (**I<sub>A</sub>**), mild (**I<sub>M</sub>**), severe (**I<sub>S</sub>**), and critical (**I<sub>C</sub>**).

## Human Synthetic Population
The human synthetic population for this ABM is assumed to have a constant size and age structure in order to reduce computational complexity. Individuals have fixed age, gender, and household assignment.
Within the ABM's daily simulation cycle, humans move between their designated home, a daytime location (either a workplace or a school), and "extracurricular" locations (other locations that people visit periodically every day as patrons).

## SARS-CoV-2 Variants of Concern
The ABM has the capacity to model the introduction and transmission of different viral variants. These variants can differ in their simulated transmissibility, pathogenicity, and immune escape capabilities (including ability to escape vaccine protection).

# Setting up the state-based counterfactual simulations
* Create a directory called "dose_data"
* From the [CDC COVID Data Tracker Vaccination Trends](https://covid.cdc.gov/covid-data-tracker/#vaccination-trends)
    - Under the "Data Table for Trends in Number of COVID-19 Vaccinations in the US" dropdown panel, there is a button to "Download Data"
    - Download data for Florida, Mississippi, and Vermont by changing the state with the "Select a Location" dropdown at the top of the screen
    - Move the downloaded files to the "dose_data" directory
* Run the following in order:
    - `bash download_cdc_covid_vax_data.sh`
    - `Rscript vax_adj_v7.R ACS_2019_pop_data.csv cdc_covid-19_vax_data.csv ./dose_data`
    - `make`
    - `./sim_test abc_covid_counterfactuals.json --process`
    - `./sim_test abc_covid_counterfactuals.json --simulate --serial 0`
---
# Setting up active vaccination counterfactual simulations
TBA

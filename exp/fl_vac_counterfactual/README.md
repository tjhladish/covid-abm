# Setting up the state-based counterfactual simulations
* Create a directory called "dose_data"
* From the [CDC COVID Data Tracker Vaccination Trends](https://covid.cdc.gov/covid-data-tracker/#vaccination-trends)
    - Under the "Data Table for Trends in Number of COVID-19 Vaccinations in the US" dropdown panel, there is a button to "Download Data"
    - Download data for Florida, Mississippi, and Vermont by changing the state with the "Select a Location" dropdown at the top of the screen
    - Move the downloaded files to the "dose_data" directory
* Run the following in order:
    - `bash download_cdc_covid_vax_data.sh`
    - `make process_state_vax_data`
    - `make`
    - `./sim_test abc_covid_counterfactuals.json --process`
    - `./sim_test abc_covid_counterfactuals.json --simulate --serial 0`
---
# Setting up active vaccination counterfactual simulations
TBA

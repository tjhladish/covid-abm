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
---
# Processing scripts
* To run simulation diagnostics/data processing:
    - In "main.cpp", set `par->dump_simulation_data` to `true`
    - After the simulation, a database named `sim_data_<serial>.sqlite` is generated
    - Run the following R scripts to analyze model performance:
        - `Rscript sim_data_processing_v1.2.R sim_data_<serial>.sqlite`
        - `Rscript model_vax_delivery_v1.2.R ACS_2019_pop_data.csv cdc_covid-19_vax_data.csv ./dose_data sim_data_<serial>.sqlite`

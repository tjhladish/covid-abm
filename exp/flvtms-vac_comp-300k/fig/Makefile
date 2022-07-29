
default: doses.png

-include local.mk

# setup definitions

DATADIR ?= ..
SIMFILES ?= $(patsubst %,covid-active-%_w-meta.sqlite,vac vac_ring-vac-alternatives)

R = $(strip Rscript $^ $(1) $@)

# reference data for generating figures
# note: these are big files

digest.rds: digest.R $(addprefix ${DATADIR}/,${SIMFILES})
	$(call R)

digest-key.rds digest-doses.rds digest-ref.rds: digest.rds

vis_support.rda: vis_support.R
	$(call R)

plotlog: ${DATADIR}/active_vac_raw_output.tgz
	mkdir $@
	tar -tf $< | grep -E "log1[0-9]{3}\." | xargs tar -xf $< -C $@

cleanplotlog:
	rm -rf plotlog

PLOTINS := digest.rds digest-key.rds digest-ref.rds vis_support.rda

# main text figures

combo-low.png: combovis.R $(PLOTINS)
	$(call R)

combo-low-mt.png: combovis_mt.R $(PLOTINS)
	$(call R)

combo-low-matched-noq.png: combovis_matched.R $(PLOTINS)
	$(call R)

combo-low-matched-isq.png: combo-low-matched-noq.png

EXTDATA := rcasedeath-florida Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Vaccination_Status\
COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries\
CDC_seroprev_long dose_data/trends_in_number_of_covid19_vaccinations_in_fl

validation.png: expanded_simvis.R validation.rds $(patsubst %,${DATADIR}/%.csv,${EXTDATA}) vis_support.rda
	$(call R)

mtfigs: combo-low.png combo-low-mt.png combo-low-matched.png

# assorted initial evaluate figures; potentially supp items?

effectiveness.png: effvis.R $(filter-out digest-ref.rds,$(PLOTINS))
	$(call R)

doses.png: dosevis.R digest-doses.rds digest-key.rds vis_support.rda
	$(call R)

averted.png: avertvis.R $(filter-out digest-ref.rds,$(PLOTINS))
	$(call R)

incidence.png: incvis.R $(PLOTINS)
	$(call R)

sifigs: doses.png effectiveness.png averted.png incidence.png

allfigs: mtfigs sifigs

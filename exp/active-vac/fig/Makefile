
default: res_fig1.png res_fig1_noq.png

-include local.mk

# setup definitions

DATADIR ?= ..
SIMFILES ?= $(patsubst %,covid-active-%.sqlite,v2.0)

R = $(strip Rscript $^ $(1) $@)

# reference data for generating figures
# note: these are big files

digest.rds: digest.R $(addprefix ${DATADIR}/,${SIMFILES})
	$(call R)

digest-key.rds digest-doses.rds digest-ref.rds: digest.rds

vis_support.rda: vis_support.R
	$(call R)

plotlog: ${DATADIR}/active_vac_logs-v2.0-FL_ref_only.tgz
	mkdir $@
	tar -xf $< --strip-components=1 -C $@

digest-ext.rds: digest-ext.R ${DATADIR}/covid-active-only.sqlite digest-key.rds digest-ref.rds
	$(call R)

cleanplotlog:
	rm -rf plotlog plotlog_ao

PLOTINS := digest.rds digest-key.rds digest-ref.rds vis_support.rda

# main text figures

active-only.png: onlyvis.R digest-ext.rds digest-ext-key.rds digest-ref.rds vis_support.rda
	$(call R)

combo-low.png: combovis.R $(PLOTINS)
	$(call R)

combo-low-mt.png: combovis_mt.R $(PLOTINS)
	$(call R)

res_fig1.png: combovis_matched.R $(PLOTINS)
	$(call R)

res_fig1_noq.png: combovis_matched_noq.R $(PLOTINS)
	$(call R)

combo-low-matched-isq.png: combo-low-matched-noq.png

EXTDATA := rcasedeath-florida Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Vaccination_Status\
COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries\
CDC_seroprev_long population-pseudo-300K state_based_counterfactual_doses

validation.rds: validation_ingest.R plotlog
	$(call R)

validation.png: expanded_simvis.R validation.rds $(patsubst %,${DATADIR}/%.csv,${EXTDATA}) vis_support.rda
	$(call R)

modelinputs.png: model_inputs.R validation.rds $(patsubst %,${DATADIR}/%.csv,${EXTDATA} reporting_dump) vis_support.rda
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

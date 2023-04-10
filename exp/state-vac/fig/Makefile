
default: fig

R = $(strip Rscript $^ $(1) $@)

-include local.make

# setup definitions

# for CABP, now overriden to point to GDRIVE fig dir
DATADIR ?= ..
RESDIR := process
OUTDIR ?= output

SIMFILES ?= $(patsubst %,covid-state-%-reserialed.sqlite,v2.0)

digest: ${RESDIR}/digest.rds
fig: $(addprefix ${OUTDIR}/,combo.png combo_alt.png state_dth_vax_comp.png)

${RESDIR}/${SIMFILES}: #| ${DATADIR}
	cp $(subst ${RESDIR},${DATADIR},$@) $@

${RESDIR}:
	mkdir $@

#${RESDIR}/outcomes.rds: ${DATADIR}/rcasedeath-florida.csv
#${RESDIR}/vaccines.rds: ${DATADIR}/Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Vaccination_Status.csv
#${RESDIR}/hospitals.rds: ${DATADIR}/COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries.csv
#${RESDIR}/seroprev.rds: ${DATADIR}/CDC_seroprev_long.csv
#${RESDIR}/FLvaccines.rds: ${DATADIR}/population-pseudo-300K.csv ${DATADIR}/state_based_counterfactual_doses.csv
#${RESDIR}/detection.rds: ${DATADIR}/reporting_dump.csv
#${RESDIR}/validation.rds: plotlog

${RESDIR}/%.rds: ${RESDIR}/%.R vis_support.rda | ${RESDIR}
	$(call R)

# reference data for generating figures
# note: these are big files
${RESDIR}/digest.rds: digest.R $(addprefix ${RESDIR}/,${SIMFILES}) | ${RESDIR}
	Rscript --min-vsize=12G $^ $@

$(patsubst %,${RESDIR}/digest-%.rds,key doses ref): ${RESDIR}/digest.rds

combo.png: ${RESDIR}/digest.rds ${RESDIR}/digest-key.rds
combo_alt.png: ${RESDIR}/digest.rds ${RESDIR}/digest-key.rds
state_dth_vax_comp.png: $(patsubst %,${RESDIR}/%.rds,acs_pop cdc_vax us_cases_deaths us_hosp)

${OUTDIR}/%.png: %.png | ${OUTDIR}
	cp $< $@

%.png: %.R vis_support.rda | ${OUTDIR}
	$(call R)

vis_support.rda: vis_support.R
	$(call R)

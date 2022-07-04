-include local.mk

DATADIR ?= ../
SIMFILES ?= abc-flvtms-vac_comp.json covid_counterfactuals_v1.0.sqlite

R = $(strip Rscript $^ $(1) $@)

${DATADIR}/digest.rds: digest.R $(addprefix ${DATADIR}/,${SIMFILES})
	$(call R)

effectiveness.png: effvis.R ${DATADIR}/digest.rds
	$(call R)


DATADIR ?= ../tmp
SIMFILES ?= abc_covid_counterfactuals.json covid_counterfactuals_v1.0.sqlite

R = $(strip Rscript $^ $(1) $@)

${DATADIR}/digest.rds: digest.R $(addprefix ${DATADIR}/,${SIMFILES})
	$(call R)

effectiveness.png: effvis.R ${DATADIR}/digest.rds
	$(call R)

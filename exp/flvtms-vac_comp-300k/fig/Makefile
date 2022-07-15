-include local.mk

DATADIR ?= ..
SIMFILES ?= abc-active-vac.json covid-active-vac_w-meta.sqlite

R = $(strip Rscript $^ $(1) $@)

digest.rds: digest.R $(addprefix ${DATADIR}/,${SIMFILES})
	$(call R)

digest-key.rds: digest.rds

vis_support.rda: vis_support.R
	$(call R)

effectiveness.png: effvis.R digest.rds digest-key.rds vis_support.rda
	$(call R)

averted.png: avertvis.R digest.rds digest-key.rds vis_support.rda
	$(call R)

incidence.png: incvis.R digest.rds digest-key.rds vis_support.rda
	$(call R)

figs: effectiveness.png averted.png incidence.png

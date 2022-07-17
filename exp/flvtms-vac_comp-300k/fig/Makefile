-include local.mk

DATADIR ?= ..
SIMFILES ?= abc-active-vac.json covid-active-vac_w-meta.sqlite

R = $(strip Rscript $^ $(1) $@)

digest.rds: digest.R $(addprefix ${DATADIR}/,${SIMFILES})
	$(call R)

digest-key.rds: digest.rds

vis_support.rda: vis_support.R
	$(call R)

PLOTINS := digest.rds digest-key.rds vis_support.rda

effectiveness.png: effvis.R $(PLOTINS)
	$(call R)

averted.png: avertvis.R $(PLOTINS)
	$(call R)

incidence.png: incvis.R $(PLOTINS)
	$(call R)

combo-low.png: combovis.R $(PLOTINS)
	$(call R)

figs: effectiveness.png averted.png incidence.png combo-low.png

-include local.mk

DATADIR ?= ..
SIMFILES ?= covid-active-vac_w-meta.sqlite covid-active-vac_ring-vac-alternatives_w-meta.sqlite

R = $(strip Rscript $^ $(1) $@)

digest.rds: digest.R $(addprefix ${DATADIR}/,${SIMFILES})
	$(call R)

digest-key.rds digest-doses.rds digest-ref.rds: digest.rds

vis_support.rda: vis_support.R
	$(call R)

PLOTINS := digest.rds digest-key.rds digest-ref.rds vis_support.rda

effectiveness.png: effvis.R $(filter-out digest-ref.rds,$(PLOTINS))
	$(call R)

averted.png: avertvis.R $(filter-out digest-ref.rds,$(PLOTINS))
	$(call R)

incidence.png: incvis.R $(PLOTINS)
	$(call R)

combo-low.png: combovis.R $(PLOTINS)
	$(call R)

combo-low-mt.png: combovis_mt.R $(PLOTINS)
	$(call R)

combo-low-matched.png: combovis_matched.R $(PLOTINS)
	$(call R)

figs: effectiveness.png averted.png incidence.png combo-low.png combo-low-mt.png combo-low-matched.png

thefig: combo-low-matched.png

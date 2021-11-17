
DEFEXT ?= png
RESDB ?= $(wildcard covid_vac_*.sqlite)
R = $(strip Rscript $^ $(1) $@)
figtar = $(addprefix $(1).,png tif)

default: $(patsubst %,%.${DEFEXT},critical_immunization_threshold variant-impact vacc-impact fl-vac-mutation-0 fl-vac-mutation-1)

rawresults.rds: injest.R ${RESDB}
	$(call R)

%eff.rds: %eff.R rawresults.rds
	$(call R)

rt.rds: rt.R rawresults.rds
	$(call R)

plottools.rda: plottools.R
	$(call R)

process: $(addsuffix .rds,vaceff vareff rt) plottools.rda

$(call figtar,critical_immunization_threshold): critical_immunization_threshold.R
	$(call R)

$(call figtar,variant-impact): variant-plots.R rawresults.rds
	$(call R)

$(call figtar,vacc-impact): vacc-plots.R rawresults.rds
	$(call R)

$(call figtar,fl-vac-mutation-%): fl-vac-plots.R vaceff.rds rt.rds plottools.rda
	$(call R,$*)

clean:
	rm *.${DEFEXT}

default: vfig

-include local.make

# setup definitions

DATADIR ?= ..
DLDIR ?= ~/Downloads
RESDIR := process
# for CABP, now overriden to point to GDRIVE fig dir
OUTDIR ?= output
SIMFILES ?= $(patsubst %,covid-active-%.sqlite,v5.1)

R = $(strip Rscript $^ $(1) $@)

vfig: $(addprefix ${OUTDIR}/,model_inputs.png validation_review.png)

${DATADIR}/${SIMFILES}:
	mv $(subst ${DATADIR},${DLDIR},$@) $@

${RESDIR}:
	mkdir $@

${RESDIR}/outcomes.rds: ${DATADIR}/rcasedeath-florida.csv
${RESDIR}/vaccines.rds: ${DATADIR}/Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Vaccination_Status.csv
${RESDIR}/hospitals.rds: ${DATADIR}/COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_State_Timeseries.csv
${RESDIR}/seroprev.rds: ${DATADIR}/CDC_seroprev_long.csv
${RESDIR}/FLvaccines.rds: ${DATADIR}/population-pseudo-300K.csv ${DATADIR}/state_based_counterfactual_doses.csv
${RESDIR}/detection.rds: ${DATADIR}/reporting_by_day.txt
${RESDIR}/vocpattern.rds: ${RESDIR}/validation.rds
${RESDIR}/vocwindows.rds: ${RESDIR}/vocpattern.rds
${RESDIR}/validation.rds: plotlog

${RESDIR}/%.rds: ${RESDIR}/%.R vis_support.rda | ${RESDIR}
	$(call R)

# reference data for generating figures
# note: these are big files
${RESDIR}/digest.rds: digest.R $(addprefix ${DATADIR}/,${SIMFILES}) | ${RESDIR}
	Rscript --min-vsize=12G $^ $@

$(patsubst %,${RESDIR}/digest-%.rds,key doses ref): ${RESDIR}/digest.rds

${RESDIR}/alt_eff.rds: alt_eff_digest.R $(addprefix ${DATADIR}/,${SIMFILES}) | ${RESDIR}
	Rscript --min-vsize=12G $^ $@

${RESDIR}/alt_q_eff.rds: alt_q_eff_digest.R $(addprefix ${DATADIR}/,${SIMFILES}) | ${RESDIR}
	Rscript --min-vsize=12G $^ $@

${RESDIR}/alt_i_eff.rds: alt_i_eff_digest.R $(addprefix ${DATADIR}/,${SIMFILES}) | ${RESDIR}
	Rscript --min-vsize=12G $^ $@

${RESDIR}/alt_eff-rt.rds: ${RESDIR}/alt_eff.rds

INPUTS := $(patsubst %,${RESDIR}/%.rds,validation outcomes vaccines hospitals seroprev FLvaccines detection vocpattern vocwindows)

allinputs: ${INPUTS}
cleaninputs:
	cd ${RESDIR} && rm *.rds

PLOTINS := $(addprefix ${RESDIR}/,digest.rds digest-key.rds digest-ref.rds)

${OUTDIR}/model_inputs.png: ${INPUTS}
${OUTDIR}/validation_review.png: ${INPUTS}
${OUTDIR}/everything.png: ${INPUTS}
${OUTDIR}/everything_bw.png: ${INPUTS}
${OUTDIR}/avertvis.png: ${PLOTINS}
${OUTDIR}/cumulative_doses.png: ${RESDIR}/digest-doses.rds ${RESDIR}/digest-key.rds
${OUTDIR}/summary.png: ${RESDIR}/alt_eff.rds ${RESDIR}/digest-key.rds
${OUTDIR}/con_summary.png: ${RESDIR}/alt_eff.rds ${RESDIR}/digest-key.rds

${OUTDIR}/%.png: %.R vis_support.rda | ${OUTDIR}
	$(call R)

${OUTDIR}/avertvis_%.png: avertvis_outcome.R vis_support.rda ${PLOTINS} | ${OUTDIR}
	$(call R)

${OUTDIR}/cum_averted_%.png: cum_avert_outcome.R vis_support.rda ${PLOTINS} | ${OUTDIR}
	$(call R)

${OUTDIR}/cum_eff_%.png: cum_eff_outcome.R vis_support.rda ${PLOTINS} | ${OUTDIR}
	$(call R)

${OUTDIR}/alt_eff_%.png: alt_eff_outcome.R vis_support.rda ${RESDIR}/alt_eff.rds \
${RESDIR}/digest-key.rds | ${OUTDIR}
	$(call R)

${OUTDIR}/alt_ave_%.png: alt_ave_outcome.R vis_support.rda ${RESDIR}/alt_eff.rds \
${RESDIR}/digest-key.rds | ${OUTDIR}
	$(call R)

${OUTDIR}/%_eff_all.png: %_eff_all.R vis_support.rda ${RESDIR}/alt_eff.rds \
${RESDIR}/digest-key.rds ${RESDIR}/vocwindows.rds | ${OUTDIR}
	$(call R)

${OUTDIR}/alt_ave_all.png: alt_ave_all.R vis_support.rda ${RESDIR}/alt_eff.rds \
${RESDIR}/digest-key.rds ${RESDIR}/vocwindows.rds | ${OUTDIR}
	$(call R)

${OUTDIR}/%_ci_all.png: %_ci_all.R vis_support.rda ${RESDIR}/alt_eff.rds \
${RESDIR}/digest-key.rds ${RESDIR}/vocwindows.rds | ${OUTDIR}
	$(call R)

${OUTDIR}/alt_mul_all.png: alt_mul_all.R vis_support.rda ${RESDIR}/alt_eff.rds \
${RESDIR}/digest-key.rds ${RESDIR}/vocwindows.rds | ${OUTDIR}
	$(call R)

${OUTDIR}/alt_inc_%.png: alt_inc_outcome.R vis_support.rda ${RESDIR}/alt_eff.rds \
${RESDIR}/digest-key.rds | ${OUTDIR}
	$(call R)

allavert: $(addprefix ${OUTDIR}/,$(patsubst %,avertvis_%.png,inf sev deaths) $(patsubst %,cum_averted_%.png,inf sev deaths) avertvis.png)
alleff: $(addprefix ${OUTDIR}/,$(patsubst %,cum_eff_%.png,inf sev deaths))
alteff: $(addprefix ${OUTDIR}/,$(patsubst %,alt_eff_%.png,inf sev deaths))
altave: $(addprefix ${OUTDIR}/,$(patsubst %,alt_ave_%.png,inf sev deaths))
altinc: $(addprefix ${OUTDIR}/,$(patsubst %,alt_inc_%.png,inf sev deaths))
allfigs: $(patsubst %,${OUTDIR}/%.png,model_input validation_review everything)

altall: $(addprefix ${OUTDIR}/,$(patsubst %,alt_%_all.png,eff ave ci inc mul) cumulative_doses.png everything.png everything_bw.png)
conall: $(addprefix ${OUTDIR}/,$(patsubst %,conditioned_%_all.png,eff ci) con_summary.png)

vis_support.rda: vis_support.R
	$(call R)

${DATADIR}/active_vac_logs-v4.0-FL_ref_only.tgz:
	mv $(subst ${DATADIR},${DLDIR},$@) $@

plotlog: | ${DATADIR}/active_vac_logs-v4.0-FL_ref_only.tgz
	mkdir $@
	tar -xf $| --strip-components=1 -C $@

${OUTDIR}/fullpage_%.png: fullpage.R vis_support.rda ${RESDIR}/digest.rds ${RESDIR}/digest-key.rds ${RESDIR}/digest-ref.rds
	$(call R)

fullpagers: $(patsubst %,${OUTDIR}/fullpage_%.png,inf_value inf_averted inf_ceff sev_value sev_averted sev_ceff deaths_value deaths_averted deaths_ceff)

digest-ext.rds: digest-ext.R ${DATADIR}/covid-active-only.sqlite digest-key.rds digest-ref.rds
	$(call R)

cleanplotlog:
	rm -rf plotlog plotlog_ao

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

validation.png: expanded_simvis.R validation.rds $(patsubst %,${DATADIR}/%.csv,${EXTDATA}) vis_support.rda
	$(call R)

mtfigs: combo-low.png combo-low-mt.png combo-low-matched.png

# assorted initial evaluate figures; potentially supp items?

effectiveness.png: effvis.R $(filter-out digest-ref.rds,$(PLOTINS))
	$(call R)

doses.png: dosevis.R ${RESDIR}/digest-doses.rds ${RESDIR}/digest-key.rds vis_support.rda
	$(call R)

incidence.png: incvis.R $(PLOTINS)
	$(call R)

sifigs: doses.png effectiveness.png averted.png incidence.png

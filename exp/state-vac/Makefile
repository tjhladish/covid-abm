-include local.mk

CPP = g++

WORKDIR ?= $(HOME)/work
CFLAGS = -std=c++2a -Wall -Wextra -Wno-deprecated-declarations --pedantic
ABCDIR = $(WORKDIR)/AbcSmc
COVDIR = $(WORKDIR)/covid-abm
GSL_PATH = $(WORKDIR)/AbcSmc/gsl_local
COVOBJ = $(COVDIR)/Person.o $(COVDIR)/Location.o $(COVDIR)/Community.o $(COVDIR)/Parameters.o $(COVDIR)/Utility.o $(COVDIR)/Vac_Campaign.o
SQLDIR = $(ABCDIR)/sqdb
DOSEDIR = $(COVDIR)/exp/active-vac/dose_data
DOSEDATA = $(DOSEDIR)/trends_in_number_of_covid19_vaccinations_in_fl.csv $(DOSEDIR)/trends_in_number_of_covid19_vaccinations_in_vt.csv $(DOSEDIR)/trends_in_number_of_covid19_vaccinations_in_ms.csv
OUTPUTDIR ?= $(COVDIR)/exp/state-vac

PATHDEFS = -DWORKDIR="\"$(WORKDIR)\"" -DOUTPUTDIR="\"$(OUTPUTDIR)\""

INCLUDE = -I$(ABCDIR) -I$(COVDIR) -I$(GSL_PATH)/include/ -I$(ABCDIR)/jsoncpp/include
#ifdef HPC_GSL_INC
#INCLUDE += -I$$HPC_GSL_INC
#endif

ABC_LIB = -L$(ABCDIR) -L$(COVDIR) -labc
GSL_LIB = -lm -L$(GSL_PATH)/lib/ -lgsl -lgslcblas -lpthread -ldl

default: libabc covid main.cpp
	$(CPP) -O2 $(CFLAGS) $(PATHDEFS) $(INCLUDE) -I$(SQLDIR) main.cpp -o sim_test $(COVOBJ) $(ABC_LIB) $(GSL_LIB)

no_db: libabc covid main-no_db.cpp
	$(CPP) -O2 $(CFLAGS) $(PATHDEFS) $(INCLUDE) main-no_db.cpp -o sim_test-no_db $(COVOBJ) $(ABC_LIB) $(GSL_LIB)

debug: libabc covid_debug main.cpp
	$(CPP) -g $(CFLAGS) $(PATHDEFS) $(INCLUDE) -I$(SQLDIR) main.cpp -o sim_test $(COVOBJ) $(ABC_LIB) $(GSL_LIB)

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

covid:
	$(MAKE) -C $(COVDIR) -f Makefile

covid_debug:
	$(MAKE) -C $(COVDIR) -f Makefile OPTI=-g

abc_sql: libabc covid main_abc.cpp
	$(CPP) -O2 $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_abc.cpp -o abc_sql $(COVOBJ) $(ABC_LIB) $(GSL_LIB)

abc_sql_debug: libabc covid main_abc.cpp
	$(CPP) -g $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_abc.cpp -o abc_sql $(COVOBJ) $(ABC_LIB) $(GSL_LIB)

process_state_vax_data: vax_adj_v7.R ACS_2019_pop_data.csv cdc_covid-19_vax_data.csv $(DOSEDATA)
	Rscript vax_adj_v7.R ACS_2019_pop_data.csv cdc_covid-19_vax_data.csv $(DOSEDIR)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(COVDIR) -f Makefile clean
	rm -f abc abc_sql

# assumes invoked from this directory...

figs:
	${MAKE} -wC fig

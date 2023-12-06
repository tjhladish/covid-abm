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
#OUTPUTDIR ?= /red/longini/tjhladish/work/covid-abm/active-vac
OUTPUTDIR ?= $(COVDIR)/exp/active-vac

PATHDEFS = -DWORKDIR="\"$(WORKDIR)\"" -DOUTPUTDIR="\"$(OUTPUTDIR)\""

INCLUDE = -I$(COVDIR)
INCLUDE += -I$(ABCDIR)/include -I$(ABCDIR)/lib -I$(ABCDIR)/lib/sqdb/include
INCLUDE += -I$(ABCDIR)/lib/CCRC32/include -I$(ABCDIR)/lib/PLS/include
INCLUDE += -I$(ABCDIR)/lib/PLS/lib/eigen -I$(ABCDIR)/lib/jsoncpp/include
#ifdef HPC_GSL_INC
#INCLUDE += -I$$HPC_GSL_INC
#endif

LIB := -L$(ABCDIR)/build -L$(ABCDIR)/build/PLS -labc -lpls -ljsoncpp -lm -lgsl -lgslcblas

default: libabc covid main.cpp
	$(CPP) -O2 $(CFLAGS) $(PATHDEFS) $(INCLUDE) -I$(SQLDIR) main.cpp -o sim_test $(COVOBJ) $(LIB)

debug: libabc covid_debug main.cpp
	$(CPP) -g $(CFLAGS) $(PATHDEFS) $(INCLUDE) -I$(SQLDIR) main.cpp -o sim_test $(COVOBJ) $(LIB)

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile abc

covid:
	$(MAKE) -C $(COVDIR) -f Makefile

covid_debug:
	$(MAKE) -C $(COVDIR) -f Makefile OPTI=-g

abc_sql: libabc covid main_abc.cpp
	$(CPP) -O2 $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_abc.cpp -o abc_sql $(COVOBJ) $(LIB)

abc_sql_debug: libabc covid main_abc.cpp
	$(CPP) -g $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_abc.cpp -o abc_sql $(COVOBJ) $(LIB)

process_state_vax_data: vax_adj_v7.R ACS_2019_pop_data.csv cdc_covid-19_vax_data.csv $(DOSEDATA)
	Rscript vax_adj_v7.R ACS_2019_pop_data.csv cdc_covid-19_vax_data.csv $(DOSEDIR)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(COVDIR) -f Makefile clean
	rm -f abc abc_sql

# assumes invoked from this directory...

figs:
	${MAKE} -wC fig

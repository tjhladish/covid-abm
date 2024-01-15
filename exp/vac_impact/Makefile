-include ../common.mk
-include local.mk

#OUTPUTDIR ?= /red/longini/tjhladish/work/covid-abm/active-vac
OUTPUTDIR ?= $(COVDIR)/exp/active-vac

PATHDEFS = -DWORKDIR="\"$(WORKDIR)\"" -DOUTPUTDIR="\"$(OUTPUTDIR)\""

#ifdef HPC_GSL_INC
#INCLUDE += -I$$HPC_GSL_INC
#endif

default: libabc covid main.cpp
	$(CPP) -O2 $(CFLAGS) $(PATHDEFS) $(INCLUDE) -I$(SQLDIR) main.cpp -o sim_test $(COVOBJ) $(LIB)

debug: libabc covid_debug main.cpp
	$(CPP) -g $(CFLAGS) $(PATHDEFS) $(INCLUDE) -I$(SQLDIR) main.cpp -o sim_test $(COVOBJ) $(LIB)

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

covid:
	$(MAKE) -C $(COVDIR) -f Makefile

covid_debug:
	$(MAKE) -C $(COVDIR) -f Makefile OPTI=-g

abc_sql: libabc covid main_abc.cpp
	$(CPP) -O2 $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_abc.cpp -o abc_sql $(COVOBJ) $(LIB)

abc_sql_debug: libabc covid main_abc.cpp
	$(CPP) -g $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_abc.cpp -o abc_sql $(COVOBJ) $(LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(COVDIR) -f Makefile clean
	rm -f abc abc_sql

# assumes invoked from this directory...

figs:
	${MAKE} -wC fig

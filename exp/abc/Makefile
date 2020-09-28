-include local.mk

CPP = g++

#CFLAGS = -g -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
CFLAGS = -O2 -std=c++11 -Wall -Wextra -Wno-deprecated-declarations --pedantic
ABCDIR = $(HOME)/work/AbcSmc
COVDIR = $(HOME)/work/covid-abm
GSL_PATH = $(HOME)/work/AbcSmc/gsl_local
COVOBJ = $(COVDIR)/Person.o $(COVDIR)/Location.o $(COVDIR)/Community.o $(COVDIR)/Parameters.o $(COVDIR)/Utility.o
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(COVDIR) -I$(GSL_PATH)/include/ -I$(ABCDIR)/jsoncpp/include
#ifdef HPC_GSL_INC
#INCLUDE += -I$$HPC_GSL_INC
#endif

ABC_LIB = -L$(ABCDIR) -L$(COVDIR) -labc -ljsoncpp -lsqdb $(ABCDIR)/sqlite3.o
GSL_LIB = -lm -L$(GSL_PATH)/lib/ -lgsl -lgslcblas -lpthread -ldl

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

covid:
	$(MAKE) -C $(COVDIR) -f Makefile

abc_sql: libabc covid main_abc.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) -I$(SQLDIR) main_abc.cpp -o abc_sql $(COVOBJ) $(ABC_LIB) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) -f Makefile clean
	$(MAKE) -C $(COVDIR) -f Makefile clean
	rm -f abc abc_sql

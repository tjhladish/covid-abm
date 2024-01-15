CPP = g++
CFLAGS = -std=c++2a -Wall -Wextra -Wno-deprecated-declarations --pedantic
WORKDIR ?= $(HOME)/work
ABCDIR = $(WORKDIR)/AbcSmc
COVDIR = $(WORKDIR)/covid-abm
COVOBJ = $(COVDIR)/Person.o $(COVDIR)/Location.o $(COVDIR)/Community.o $(COVDIR)/Parameters.o $(COVDIR)/Utility.o $(COVDIR)/Vac_Campaign.o
SQLDIR = $(ABCDIR)/sqdb

INCLUDE = -I$(ABCDIR) -I$(DENDIR) -I$(IMMDIR) -I$(ABCDIR)/jsoncpp/include

INCLUDE = -I$(COVDIR)
INCLUDE += -I$(ABCDIR)/include -I$(ABCDIR)/lib -I$(ABCDIR)/lib/sqdb/include
INCLUDE += -I$(ABCDIR)/lib/CCRC32/include -I$(ABCDIR)/lib/PLS/include
INCLUDE += -I$(ABCDIR)/lib/PLS/lib/eigen -I$(ABCDIR)/lib/jsoncpp/include

#ifdef TACC_GSL_INC
#INCLUDE += -I$$TACC_GSL_INC
#endif

ifdef HPC_GSL_INC
INCLUDE += -I$$HPC_GSL_INC
endif

LIB := -L$(ABCDIR)/build -L$(ABCDIR)/build/PLS -labc -lpls -ljsoncpp -lm -lgsl -lgslcblas -lpthread -ldl

ifdef HPC_GSL_LIB
LIB += -L$$HPC_GSL_LIB
endif


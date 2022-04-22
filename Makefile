SHELL=/bin/bash
G++VER := $(shell command -v g++-4.9)

ifndef G++VER
CPP:=g++
else
CPP:=g++-4.9
endif

GSL_PATH = $(HOME)/work/AbcSmc/gsl_local

CFLAGS   	= -Wall -Wextra -pedantic -std=c++17
OPTI     	?= -O2
LDFLAGS	 	= -L$(GSL_PATH)/lib/ # $(HPC_GSL_LIB) $(TACC_GSL_LIB)
INCLUDES 	= -I$(GSL_PATH)/include # $(HPC_GSL_INC) $(TACC_GSL_INC)
LIBS     	= -lm -lgsl -lgslcblas
DEFINES  	= -DVERBOSE

default: Makefile Vac_Campaign.o Person.o Location.o Community.o Parameters.o Utility.o

%.o: %.cpp Community.h Location.h Utility.h Parameters.h Person.h Vac_Campaign.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c $<

date_test: date_test.cpp Parameters.o Utility.o Date.h
	$(CPP) $(CFLAGS) $(OPTI) date_test.cpp -o date_test /home/tjhladish/work/covid-abm/Parameters.o /home/tjhladish/work/covid-abm/Utility.o -lgsl

clean:
	rm -f *.o model *~

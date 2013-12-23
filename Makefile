# Makefile for RNAmfpt

CCFLAGS           = -c -std=c99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q -I $(H)
LDFLAGS           = -lm -lgomp -llapack -L/usr/local/include -llapacke -lgslcblas -lgsl -o
BINDIR            = ~/bin # Change this to the BINDIR
LIBDIR            = ~/lib # Change this to the LIBDIR
CC                = gcc
GCC_VERSION      := $(shell expr `$(CC) -dumpversion`)
CC_MAJ_VERSION   := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 1` \* 10000)
CC_MIN_VERSION   := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 2` \* 100)
CC_PATCH_VERSION := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 3`)
GCC_NUM_VERSION  := $(shell expr $(CC_MAJ_VERSION) \+ $(CC_MIN_VERSION) \+ $(CC_PATCH_VERSION))
GCC_GTEQ_4.6.0   := $(shell expr $(GCC_NUM_VERSION) \>= 40600)
GCC_GTEQ_4.9.0   := $(shell expr $(GCC_NUM_VERSION) \>= 40900)
LIB              := ../../lib
H                := ../../h

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	CCFLAGS += -Ofast -march=native
else
	CCFLAGS += -O3
endif

ifeq "$(GCC_GTEQ_4.9.0)" "1"
	CCFLAGS += -fdiagnostics-color=always
endif

RNAmfpt.out: rna_mfpt.o mfpt_params.o mfpt_energy_grid.o mfpt_parser.o
	$(CC) rna_mfpt.o mfpt_params.o mfpt_energy_grid.o $(LDFLAGS) RNAmfpt.out
	ar cr $(LIB)/libmfpt.a mfpt_parser.o mfpt_energy_grid.o mfpt_params.o
	
rna_mfpt.o: rna_mfpt.c $(H)/mfpt_parser.h $(H)/mfpt_energy_grid.h $(H)/constants.h $(H)/mfpt_params.h
	$(CC) $(CCFLAGS) rna_mfpt.c

mfpt_parser.o: mfpt_parser.c $(H)/mfpt_parser.h
	$(CC) $(CCFLAGS) mfpt_parser.c
	
mfpt_params.o: mfpt_params.c $(H)/mfpt_params.h
	$(CC) $(CCFLAGS) mfpt_params.c
	
mfpt_energy_grid.o: mfpt_energy_grid.c $(H)/mfpt_energy_grid.h $(H)/constants.h $(H)/mfpt_params.h
	$(CC) $(CCFLAGS) mfpt_energy_grid.c

clean:
	rm -f *.o $(LIB)/libmfpt.a RNAmfpt.out

install: RNAmfpt.out
	cp RNAmfpt $(BINDIR)/RNAmfpt
	cp $(LIB)/libmfpt.a $(LIBDIR)
	
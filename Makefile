# Makefile for RNAmfpt

CCFLAGS          = -c -std=gnu99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q -I $(HEADER) -I $(SHARED_HEADER)
LDFLAGS          = -L . -L $(LIB)/ -L /usr/local/include -lm -lgomp -llapack -lgslcblas -lgsl -o
BINDIR           = ~/bin
LIBDIR           = ~/lib
CC               = gcc
LIB              = ../../lib
SHARED_HEADER    = ../../h
HEADER           = h
CODE             = c
GCC_VERSION      = $(shell expr `$(CC) -dumpversion`)
CC_MAJ_VERSION   = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 1` \* 10000)
CC_MIN_VERSION   = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 2` \* 100)
CC_PATCH_VERSION = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 3`)
GCC_NUM_VERSION  = $(shell expr $(CC_MAJ_VERSION) \+ $(CC_MIN_VERSION) \+ $(CC_PATCH_VERSION))
GCC_GTEQ_4.6.0   = $(shell expr $(GCC_NUM_VERSION) \>= 40600)
GCC_GTEQ_4.9.0   = $(shell expr $(GCC_NUM_VERSION) \>= 40900)

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	CCFLAGS += -Ofast -march=native
else
	CCFLAGS += -O3
endif

ifeq "$(GCC_GTEQ_4.9.0)" "1"
	CCFLAGS += -fdiagnostics-color=always
endif

RNAmfpt.out: $(CODE)/rna_mfpt.o $(CODE)/mfpt_params.o $(CODE)/mfpt_energy_grid.o $(CODE)/mfpt_parser.o $(CODE)/mfpt_initializers.o
	$(CC) $(CODE)/rna_mfpt.o $(CODE)/mfpt_parser.o $(CODE)/mfpt_params.o $(CODE)/mfpt_energy_grid.o $(CODE)/mfpt_initializers.o $(LDFLAGS) RNAmfpt.out
	ar cr $(LIB)/libmfpt.a $(CODE)/mfpt_energy_grid.o $(CODE)/mfpt_params.o $(CODE)/mfpt_initializers.o
	
$(CODE)/rna_mfpt.o: $(CODE)/rna_mfpt.c $(HEADER)/parser.h $(HEADER)/energy_grid.h $(HEADER)/constants.h $(HEADER)/params.h $(HEADER)/initializers.h
	$(CC) $(CCFLAGS) $(CODE)/rna_mfpt.c -o $(CODE)/rna_mfpt.o

$(CODE)/mfpt_initializers.o: $(CODE)/mfpt_initializers.c $(HEADER)/initializers.h
	$(CC) $(CCFLAGS) $(CODE)/mfpt_initializers.c -o $(CODE)/mfpt_initializers.o

$(CODE)/mfpt_parser.o: $(CODE)/mfpt_parser.c $(HEADER)/parser.h
	$(CC) $(CCFLAGS) $(CODE)/mfpt_parser.c -o $(CODE)/mfpt_parser.o
	
$(CODE)/mfpt_params.o: $(CODE)/mfpt_params.c $(HEADER)/params.h
	$(CC) $(CCFLAGS) $(CODE)/mfpt_params.c -o $(CODE)/mfpt_params.o
	
$(CODE)/mfpt_energy_grid.o: $(CODE)/mfpt_energy_grid.c $(HEADER)/energy_grid.h $(HEADER)/constants.h $(HEADER)/params.h $(HEADER)/initializers.h
	$(CC) $(CCFLAGS) $(CODE)/mfpt_energy_grid.c -o $(CODE)/mfpt_energy_grid.o

clean:
	rm -f $(CODE)/*.o $(LIB)/libmfpt.a RNAmfpt.out

install: RNAmfpt.out
	cp RNAmfpt $(BINDIR)/RNAmfpt
	cp $(LIB)/libmfpt.a $(LIBDIR)
	

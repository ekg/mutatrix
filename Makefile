VCFLIB_ROOT=vcflib
FASTAHACK_ROOT=fastahack
INCLUDES = -lm -L. -L$(VCFLIB_ROOT)/tabixpp/ -ltabix -std=c++0x -lz 

TABIX_OBJECTS=$(VCFLIB_ROOT)/tabixpp/tabix.o $(VCFLIB_ROOT)/tabixpp/bgzf.o
SMITHWATERMAN_OBJECTS=$(VCFLIB_ROOT)/smithwaterman/SmithWatermanGotoh.o

all: mutatrix

clean:
	rm -f mutatrix
	cd $(VCFLIB_ROOT) && $(MAKE) clean
	cd $(FASTAHACK_ROOT) && $(MAKE) clean

.PHONY: all clean

# builds vcflib
$(VCFLIB_ROOT)/Variant.o:
	cd $(VCFLIB_ROOT) && $(MAKE)

# builds fastahack
$(FASTAHACK_ROOT)/Fasta.o:
	cd $(FASTAHACK_ROOT) && $(MAKE)

mutatrix: mutatrix.cpp split.cpp Repeats.cpp $(VCFLIB_ROOT)/Variant.o $(FASTAHACK_ROOT)/Fasta.o
	g++ mutatrix.cpp split.cpp Repeats.cpp $(VCFLIB_ROOT)/Variant.o \
		$(VCFLIB_ROOT)/smithwaterman/IndelAllele.o \
		$(VCFLIB_ROOT)/smithwaterman/LeftAlign.o \
		$(VCFLIB_ROOT)/smithwaterman/disorder.c \
		$(TABIX_OBJECTS) $(SMITHWATERMAN_OBJECTS) \
		$(FASTAHACK_ROOT)/Fasta.o -o mutatrix $(INCLUDES)

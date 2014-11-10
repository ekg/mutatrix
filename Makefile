VCFLIB_ROOT=vcflib
FASTAHACK_ROOT=fastahack
INCLUDES = -I$(VCFLIB_ROOT)/src/ -I$(VCFLIB_ROOT)/ -L. -L$(VCFLIB_ROOT)/tabixpp/ -ltabix -L$(VCFLIB_ROOT)/ -lvcflib -lm -lz -std=c++0x

TABIX_OBJECTS=$(VCFLIB_ROOT)/tabixpp/tabix.o $(VCFLIB_ROOT)/tabixpp/bgzf.o
SMITHWATERMAN_OBJECTS=$(VCFLIB_ROOT)/smithwaterman/SmithWatermanGotoh.o

all: mutatrix

clean:
	rm -f mutatrix
	cd $(VCFLIB_ROOT) && $(MAKE) clean
	cd $(FASTAHACK_ROOT) && $(MAKE) clean

.PHONY: all clean

# builds vcflib
$(VCFLIB_ROOT)/libvcflib.a:
	cd $(VCFLIB_ROOT) && $(MAKE) libvcflib.a

# builds fastahack
$(FASTAHACK_ROOT)/Fasta.o:
	cd $(FASTAHACK_ROOT) && $(MAKE)

mutatrix: mutatrix.cpp split.cpp Repeats.cpp $(VCFLIB_ROOT)/libvcflib.a $(FASTAHACK_ROOT)/Fasta.o
	g++ mutatrix.cpp \
		$(FASTAHACK_ROOT)/Fasta.o -o mutatrix $(INCLUDES)

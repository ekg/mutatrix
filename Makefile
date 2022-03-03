VCFLIB_ROOT=vcflib
FASTAHACK_ROOT=fastahack
INCLUDES = -I$(VCFLIB_ROOT)/src/ -I$(VCFLIB_ROOT)/filevercmp -I$(VCFLIB_ROOT)/fastahack -I$(VCFLIB_ROOT)/ -I$(VCFLIB_ROOT)/multichoose -I$(VCFLIB_ROOT)/smithwaterman -L. -L$(VCFLIB_ROOT)/tabixpp/ -ltabixpp -L$(VCFLIB_ROOT)/build -lvcflib -lm -lz -std=c++0x

TABIX_OBJECTS=$(VCFLIB_ROOT)/tabixpp/tabix.o $(VCFLIB_ROOT)/tabixpp/bgzf.o
SMITHWATERMAN_OBJECTS=$(VCFLIB_ROOT)/smithwaterman/SmithWatermanGotoh.o

all: mutatrix

clean:
	rm -f mutatrix
	cd $(VCFLIB_ROOT) && $(MAKE) clean
	cd $(FASTAHACK_ROOT) && $(MAKE) clean

.PHONY: all clean

# builds vcflib
$(VCFLIB_ROOT)/build/libvcflib.a:
	cd $(VCFLIB_ROOT) && mkdir -p build && cd build && cmake .. && $(MAKE)

# builds fastahack
$(FASTAHACK_ROOT)/Fasta.o:
	cd $(FASTAHACK_ROOT) && $(MAKE)

mutatrix: mutatrix.cpp split.cpp Repeats.cpp $(VCFLIB_ROOT)/build/libvcflib.a $(FASTAHACK_ROOT)/Fasta.o
	g++ mutatrix.cpp \
		$(FASTAHACK_ROOT)/Fasta.o -o mutatrix $(INCLUDES)

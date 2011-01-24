all: mutatrix

clean:
	rm mutatrix

.PHONY: all clean

mutatrix: mutatrix.cpp Fasta.cpp split.cpp
	g++ mutatrix.cpp Fasta.cpp split.cpp -o mutatrix


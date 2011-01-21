#include "Fasta.h"
#include <getopt.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <time.h>
#include "mt19937ar.h"
#include <math.h>
#include "convert.h"

void printSummary() {
    cerr 
         << endl
         << "program: mutatrix (population genome simulator)" << endl
         << endl
         << "usage: mutatrix [options] reference.fa >mutants.vcf" << endl
         << endl
         << "options:" << endl 
         << "    -r, --rate             the relative rate of mutation per bp per chrom (default 0.001)" << endl
         << "    -a, --afs-alpha        the allele frequency spectrum distribution alpha parameter (zeta(n), default 1.7)" << endl
         << "    -z, --indel-alpha      the alpha parameter of the indel length frequency distribution (zeta(l), default 1.7)" << endl
         << "                           indels of length N have probability zeta(N)" << endl
         << "    -s, --indel-snp-ratio  ratio between 1bp indels and snps (default 0.2)" << endl
         << "    -M, --indel-max        maximum indel length (default 1000)" << endl
         << "    -m, --mnp-ratio        the geometric scaling probability for 2bp multi-nucleotide-polymorphisms relative to SNPs," << endl
         << "                           2bp MNPs relative to 3pb MNPs, etc. (default 0.01)" << endl
         << "    -p, --ploidy           ploidy of the population (default 1)" << endl
         << "    -n, --population-size  number of individuals in the population" << endl
         << "    -P, --output-prefix    prefix output fasta files with this" << endl
         << endl
         << "Generates a simulated population with no linkage, but a zeta-distributed allele frequency spectrum." << endl
         << "Writes a set of files PREFIX_SEQUENCE_INDIVIDUAL_COPY.fa for each fasta sequence in the provided" << endl
         << "reference, sample, and simulated copy of the genome.  A VCF file is generated on stdout describing" << endl
         << "the reference-relative variation of each sample." << endl
         << endl
         << "The allele frequency spectrum and indel length distribution are zeta-distributed.  The MNP length" << endl
         << "frequency spectrum is geometrically distributed." << endl
         << endl
         << "author: Erik Garrison <erik.garrison@bc.edu>" << endl
         << endl;
}

void writeFasta(ostream& out, string& seqname, string& sequence, int linewidth = 80) {
    out << ">" << seqname << endl;
    long int pos = 0;
    while (pos < sequence.length()) {
        out << sequence.substr(pos, linewidth) << endl;
        pos += linewidth;
    }
}


string dateStr(void) {

    time_t rawtime;
    struct tm* timeinfo;
    char buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 80, "%Y%m%d", timeinfo);

    return string(buffer);

}

int expovariate(double lambda) {
    return -log(genrand_real1()) / lambda;
}

double zetarandom(double alpha) {

    double u, v;
    double X, T;
    double test;
    double b = pow(2.0, alpha - 1.0);

    do {
        u = genrand_res53();
        v = genrand_res53();
        X = floor (pow (u, -1.0 / (alpha - 1.0)));
        T = pow (1.0 + 1.0 / X, alpha - 1.0);
        test = v * X * (T - 1.0) / (b - 1.0);
    } while ( test > (T / b) );

    return X;

}

int main (int argc, char** argv) {

    float mutation_rate = 0.001;
    float het_rate = 0.5;
    float afs_alpha = 1.7;
    float indel_alpha = 1.7;
    float indel_snp_ratio = 0.2;
    float mnp_ratio = 0.01;
    int indel_max = 1000;
    int ploidy = 1;
    int population_size = 1;
    string fastaFileName;

    string command_line = argv[0];
    for (int i = 1; i < argc; ++i) {
        command_line += " ";
        command_line += argv[i];
    }

    int c;

    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            //{"brief",   no_argument,       &verbose_flag, 0},
            {"help", no_argument, 0, 'h'},
            {"rate",  required_argument, 0, 'r'},
            {"afs-alpha",  required_argument, 0, 'a'},
            {"indel-alpha", required_argument, 0, 'z'},
            {"indel-snp-ratio", required_argument, 0, 's'},
            {"indel-max", required_argument, 0, 'M'},
            {"mnp-ratio", required_argument, 0, 'M'},
            {"ploidy", required_argument, 0, 'p'},
            {"population-size", required_argument, 0, 'n'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hr:a:z:s:p:n:M:m:", long_options, &option_index);

      /* Detect the end of the options. */
          if (c == -1)
            break;
 
          switch (c)
            {
            case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
              break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
              printf (" with arg %s", optarg);
            printf ("\n");
            break;

          case 'r':
            if (!convert(optarg, mutation_rate)) {
                cerr << "could not read -r, --rate" << endl;
                exit(1);
            }
            break;

          case 'a':
            if (!convert(optarg, afs_alpha)) {
                cerr << "could not read -a, --afs-alpha" << endl;
                exit(1);
            }
            break;
 
          case 'z':
            if (!convert(optarg, indel_alpha)) {
                cerr << "could not read -z, --indel-alpha" << endl;
                exit(1);
            }
            break;

          case 's':
            if (!convert(optarg, indel_snp_ratio)) {
                cerr << "could not read -s, --indel-snp-ratio" << endl;
                exit(1);
            }
            break;

          case 'M':
            if (!convert(optarg, indel_max)) {
                cerr << "could not read -M, --indel-max" << endl;
                exit(1);
            }
            break;
 
          case 'm':
            if (!convert(optarg, mnp_ratio)) {
                cerr << "could not read -m, --mnp-ratio" << endl;
                exit(1);
            }
            break;
 
          case 'p':
            if (!convert(optarg, ploidy)) {
                cerr << "could not read -p, --ploidy" << endl;
                exit(1);
            }
            break;
 
          case 'n':
            if (!convert(optarg, population_size)) {
                cerr << "could not read -n, --population-size" << endl;
                exit(1);
            }
            break;

          case 'h':
            printSummary();
            exit(0);
            break;
 
          case '?':
            /* getopt_long already printed an error message. */
            printSummary();
            exit(1);
            break;
 
          default:
            abort ();
          }
      }

    /* Print any remaining command line arguments (not options). */
    if (optind < argc) {
        //cerr << "fasta file: " << argv[optind] << endl;
        fastaFileName = argv[optind];
    } else {
        cerr << "please specify a fasta file" << endl;
        printSummary();
        exit(1);
    }

    int seed = time(NULL);
    init_genrand(seed); // seed mt with current time

    string seqname;
    string sequence;  // holds sequence so we can process it

    FastaReference fr(fastaFileName);

    string bases = "ATGC";

    // write the VCF header
    cout 
        << "##fileformat=VCFv4.0" << endl
        << "##fileDate=" << dateStr() << endl
        << "##source=mutate population genome simulator" << endl
        << "##seed=" << seed << endl
        << "##reference=" << fastaFileName << endl
        << "##phasing=true" << endl
        << "##commandline=" << command_line << endl
        << "##INFO=<ID=SNP,Number=0,Type=Flag,Description=\"SNP allele\">" << endl
        << "##INFO=<ID=MNP,Number=0,Type=Integer,Description=\"Length of MNP allele, if present\">" << endl
        << "##INFO=<ID=INS,Number=1,Type=Integer,Description=\"Length of insertion allele, if present\">" << endl
        << "##INFO=<ID=DEL,Number=1,Type=Integer,Description=\"Length of deletion allele, if present\">" << endl
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (int i = 0; i < population_size; ++i)
        cout << "\t" << i;
    cout << endl;

    int copies = ploidy * population_size;

    map<string, vector<string> > sequencesByRefseq;

    for (FastaIndex::iterator s = fr.index->begin(); s != fr.index->end(); ++s) {

        FastaIndexEntry& indexEntry = s->second;
        seqname = indexEntry.name;
        sequence = fr.getSequence(s->first);

        vector<string>& sequences = sequencesByRefseq[seqname];
        sequences.resize(copies);

        vector<string> genotypes;
        genotypes.resize(population_size);

        long int pos = 0;
        string ref = "";

        while (pos < sequence.size()) {

            pos += ref.size(); // step by the size of the last event
            ref = sequence.substr(pos, 1); // by default, ref is just the current base

            // flip a coin
            double mut = genrand_real1();
            // scale it into the amount of sequence at this site
            mut /= copies;
            // have we drawn a mutation event?
            if (mut < mutation_rate) {

                //cerr << seqname << "\t" << pos << "\t";

                string alt;

                int len = 1;

                bool insertion = false;
                bool deletion = false;
                bool mnp = false;
                bool snp = (genrand_real1() > indel_snp_ratio);
                // make an alternate allele
                if (snp) {

                    alt = ref;
                    while (alt == ref) {
                        alt = string(1, bases.at(genrand_int32() % 4));
                    }

                    mnp = (genrand_real1() < mnp_ratio);
                    if (mnp) {
                        snp = false;
                        int i = 1;
                        do {
                            ref += sequence.substr(pos + i, 1);
                            alt += sequence.substr(pos + i, 1);
                            ++i;
                            while (alt.at(alt.size() - 1) == ref.at(ref.size() - 1)) {
                                alt.at(alt.size() - 1) = bases.at(genrand_int32() % 4);
                            }
                        } while (genrand_real1() < mnp_ratio);
                        len = alt.size();
                    }
                // indel case
                } else {
                    // how many bp?
                    len = (int) floor(zetarandom(indel_alpha));
                    // guard against out-of-sequence indels
                    if (pos + len < sequence.size() && len <= indel_max) {
                        if (genrand_int32() % 2 == 0) {
                        // deletion
                            deletion = true;
                            ref = sequence.substr(pos, 1 + len);
                            alt = sequence.substr(pos, 1);
                        } else {
                        // insertion?
                            insertion = true;
                            // TODO ... vntr?
                            // insert some random de novo bases
                            while (alt.length() < len) {
                                alt += string(1, bases.at(genrand_int32() % 4));
                            }
                        }
                    }
                }

                string genotype;

                vector<bool> alts;
                
                // AFS simulation
                int allele_freq = min((double) copies, zetarandom(afs_alpha));

                {
                    int i = 0;
                    for (; i < allele_freq; ++i) {
                        alts.push_back(true);
                    }

                    for (; i < copies; ++i) {
                        alts.push_back(false);
                    }
                }

                // shuffle the minor alleles around the population
                random_shuffle(alts.begin(), alts.end());

                // and encode the genotypes for each individual
                for (int j = 0; j < population_size; ++j) {
                    string genotype;
                    for (int i = 0; i < ploidy; ++i) {
                        int l = (j * ploidy) + i;
                        // for simplicity, only allow reference-relative events
                        // when we are outside of the last deletion
                        if (alts.at(l)) {
                            genotype += "1|";
                            sequences.at(l) += alt;
                        } else {
                            genotype += "0|";
                            sequences.at(l) += ref;
                        }
                    }
                    genotype = genotype.substr(0, genotype.size() - 1);
                    // and record it
                    genotypes.at(j) = genotype;
                }

                // and write a line of VCF output
                // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
                cout << seqname
                     << "\t" << pos
                     << "\t" << "."
                     << "\t" << ref
                     << "\t" << alt
                     << "\t" << 99
                     << "\t" << "."
                     << "\t" << "NS=" << population_size << ";AC=" << allele_freq;
                if (snp) {
                    cout << ";SNP";
                } else if (mnp) {
                    cout << ";MNP=" << len;
                } else if (insertion) {
                    cout << ";INS=" << len;
                } else if (deletion) {
                    cout << ";DEL=" << len;
                }
                cout << "\t" << "GT";

                for (vector<string>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
                    cout << "\t" << *g;
                }
                cout << endl;

            }
            // no mutation
            else {
                for (int i = 0; i < copies; ++i) {
                    sequences.at(i) += ref;
                }
            }
        }
    }

    // print sequences to files
    // named by sample number and sequence copy
    for (int i = 0; i < population_size; ++i) {
        for (int j = 0; j < ploidy; ++j) {
            int l = (j * ploidy) + i;
            stringstream num;
            num << i << "__" << j << "__";
            for (map<string, vector<string> >::iterator s = sequencesByRefseq.begin(); s != sequencesByRefseq.end(); ++s) {
                string fullname = num.str() + s->first;
                //writeFasta(cout, fullname, s->second.at(l));
            }
        }
    }

    return 0;

}

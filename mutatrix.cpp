#include "fastahack/Fasta.h"
#include <getopt.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <time.h>
#include "mt19937ar.h"
#include <cmath>
#include <assert.h>
#include "convert.h"
#include <iomanip>
#include "Repeats.h"
#include "vcflib/Variant.h"
#include "split.h"
#include "join.h"

class Allele {
public:
    friend bool operator==(const Allele&, const Allele&);
    friend bool operator<(const Allele&, const Allele&);
    friend ostream& operator<<(ostream&, const Allele&);
    string ref;
    string alt;
    string type;
    Allele(const string& r, const string& a, const string& t = "")
        : ref(r), alt(a), type(t) { }
};

bool operator==(const Allele& a, const Allele& b) {
    return a.ref == b.ref && a.alt == b.alt;
}

bool operator<(const Allele& a, const Allele& b) {
    return convert(a) < convert(b);
}

ostream& operator<<(ostream& out, const Allele& allele) {
    out << allele.ref << "/" << allele.alt;
    return out;
}

class SampleFastaFile {

public:

    ofstream fastafile;
    long int pos;
    string linebuffer;
    string filename;
    string seqname;
    int linewidth;

    void write(string sequence) {
        linebuffer += sequence;
        while (linebuffer.length() > linewidth) {
            fastafile << linebuffer.substr(0, linewidth) << endl;
            linebuffer = linebuffer.substr(linewidth);
        }
    }

    //SampleFastaFile(void) { }

    // TODO add filname and seqname to this
    SampleFastaFile(string& m_filename, string& m_seqname, int m_linewidth = 80)
        : filename(m_filename)
        , seqname(m_seqname)
        , pos(0)
        , linewidth(m_linewidth)
    {
        fastafile.open(filename.c_str());
        if (!fastafile.is_open()) {
            cerr << "could not open " << filename << " for writing, exiting" << endl;
            exit(1);
        }
        fastafile << ">" << seqname << endl;
    }

    ~SampleFastaFile(void) {
        write(""); // flush
        fastafile << linebuffer << endl;
        fastafile.close();
    }

};

/*
map<string, int> repeat_counts(
        unsigned int position,
        string& sequence,
        int maxsize) {
    map<string, int> counts;
    for (int i = 1; i <= maxsize; ++i) {
        // subseq here i bases
        string seq = sequence.substr(position, i);
        // go left.
        int j = position - i;
        int left = 0;
        while (j - i >= 0 && seq == sequence.substr(j, i)) {
            j -= i;
            ++left;
        }
        // go right.
        j = position + i;
        int right = 0;
        while (j + i < sequence.size() && seq == sequence.substr(j, i)) {
            j += i;
            ++right;
        }
        // if we went left and right a non-zero number of times, 
        if (right > 0 || left > 0) {
            counts[seq] = right + left + 1;
        }
    }

    // filter out redundant repeat information
    if (counts.size() > 1) {
        map<string, int> filteredcounts;
        map<string, int>::iterator c = counts.begin();
        string prev = c->first;
        filteredcounts[prev] = c->second;  // shortest sequence
        ++c;
        for (; c != counts.end(); ++c) {
            int i = 0;
            string seq = c->first;
            while (i + prev.length() <= seq.length() && seq.substr(i, prev.length()) == prev) {
                i += prev.length();
            }
            if (i < seq.length()) {
                filteredcounts[seq] = c->second;
                prev = seq;
            }
        }
        return filteredcounts;
    } else {
        return counts;
    }
}
*/

void printSummary() {
    cerr 
         << endl
         << "program: mutatrix (population genome simulator)" << endl
         << endl
         << "usage: mutatrix [options] reference.fa >mutants.vcf" << endl
         << endl
         << "options:" << endl 
         << "    -s, --snp-rate          the relative rate of point mutation per bp per chrom (default 0.00001)" << endl
         << "    -M, --mnp-ratio         the geometric scaling probability for 2bp multi-nucleotide-polymorphisms relative to SNPs," << endl
         << "                            2bp MNPs relative to 3pb MNPs, etc. (default 0.01)" << endl
         << "    -i, --indel-rate        the rate of non-repeat indel mutations per bp per chrom (default 0.000001)" << endl
         << "    -X, --indel-max         maximum indel length (default 1000)" << endl
         << "    -z, --indel-alpha       the alpha parameter of the indel length frequency distribution (zeta(l), default 1.1)" << endl
         << "                            indels of length N have probability zeta(N)" << endl
         << "    -q, --repeat-max-size   maximum size of exect repeat unit in the genome to detect (default 20)" << endl
         << "    -m, --microsat-rate     the rate of microsatellite mutation at microsatellite sites (default 0.000005)" << endl
         << "    -t, --microsat-afs-alpha    alpha parameter for microsatellite allele count (default 1.7)" << endl
         << "    -j, --microsat-len-alpha    alpha parameter for microsatellite mutation length (default 1.7)" << endl
         << "    -m, --microsat-min-len  the minimum number of bases in a repeat to consider it a microsatellite (default 1)" << endl
         << "    -a, --afs-alpha         the allele frequency spectrum distribution scaling parameter (1/i * alpha, default 1.0)" << endl
         << "    -p, --ploidy            ploidy of the population (default 1)" << endl
         << "    -n, --population-size   number of individuals in the population" << endl
         << "    -P, --file-prefix       prefix output fasta files with this" << endl
         << "    -S, --sample-prefix     prefix sample names (numbers by default) with this" << endl
         << "    -g, --random-seed       provide the seed for pseudorandom generation (default, seconds since 1970)" << endl
         << "    -d, --dry-run           don't write any fasta output files, just print VCF output" << endl
         << endl
         << "Generates a simulated population with no linkage, but allele frequency spectrum drawn from 1/n," << endl
         << "where n is the minor allele frequency." << endl
         << endl
         << "Writes a set of files of the form prefix_sequence_individual_copy.fa for each fasta sequence in" << endl
         << "the provided reference, sample, and simulated copy of the genome.  A VCF file is generated on stdout"<< endl
         << "describing the reference-relative variation of each sample." << endl
         << endl
         << "The indel length distribution is zeta-distributed.  The MNP length frequency spectrum is" << endl
         << "geometrically distributed." << endl
         << endl
         << "At runtime the genome is analyzed for repeats up to a certain number of bp (default 20)." << endl
         << "If repeats are found, mutations are generated from them using the microsatellite paramaters." << endl
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

long double microsatelliteInsProb(int count) {
    return min(1.0, 3.1 * pow(10, -6) * exp(0.2 * count));
}

long double microsatelliteDelProb(int count) {
    return min(1.0, 4.0 * pow(10, -7) * exp(0.302 * count));
}

int expovariate(double lambda) {
    return -log(genrand_real1()) / lambda;
}

// generates a random allele frequency in 1/i scaled by alpha
// bounded by the number of copies at the locus
int random_allele_frequency(int copies, double alpha) {
    return min((int) floor(1 / genrand_real1() * alpha), copies);
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

    double snp_mutation_rate = 0.000025;
    double indel_mutation_rate = 0.000005;
    double het_rate = 0.5;
    double afs_alpha = 1;
    double indel_alpha = 3;
    double microsatellite_afs_alpha = 1;
    double microsatellite_len_alpha = 1.7;
    double microsatellite_mutation_rate = 0.000005;
    double mnp_ratio = 0.01;
    int microsatellite_min_length = 1;
    int indel_max = 1000;
    int ploidy = 1;
    int population_size = 1;
    int sample_id_max_digits = 1;
    int seed = time(NULL);
    string fastaFileName;
    string file_prefix = "";
    string sample_prefix = "";
    bool dry_run = false;
    int repeat_size_max = 20;

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
            {"snp-rate",  required_argument, 0, 's'},
            {"mnp-ratio", required_argument, 0, 'M'},
            {"indel-rate",  required_argument, 0, 'i'},
            {"indel-alpha", required_argument, 0, 'z'},
            {"indel-max", required_argument, 0, 'X'},
            {"repeat-size-max", required_argument, 0, 'q'},
            {"microsat-rate",  required_argument, 0, 'm'},
            {"microsat-afs-alpha", required_argument, 0, 't'},
            {"microsat-len-alpha", required_argument, 0, 'j'},
            {"microsat-min-len", required_argument, 0, 'l'},
            {"afs-alpha",  required_argument, 0, 'a'},
            {"ploidy", required_argument, 0, 'p'},
            {"population-size", required_argument, 0, 'n'},
            {"file-prefix", required_argument, 0, 'P'},
            {"sample-prefix", required_argument, 0, 'S'},
            {"random-seed", required_argument, 0, 'g'},
            {"dry-run", no_argument, 0, 'd'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "hda:z:s:i:q:p:n:M:X:t:m:P:S:g:l:j:", long_options, &option_index);

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

          case 'd':
            dry_run = true;
            break;

          case 'q':
            if (!convert(optarg, repeat_size_max)) {
                cerr << "could not read -q, --repeat-size-max" << endl;
                exit(1);
            }
            break;

          case 's':
            if (!convert(optarg, snp_mutation_rate)) {
                cerr << "could not read -s, --snp-rate" << endl;
                exit(1);
            }
            break;

          case 'i':
            if (!convert(optarg, indel_mutation_rate)) {
                cerr << "could not read -i, --indel-rate" << endl;
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

          case 'X':
            if (!convert(optarg, indel_max)) {
                cerr << "could not read -M, --indel-max" << endl;
                exit(1);
            }
            break;
 
          case 'M':
            if (!convert(optarg, mnp_ratio)) {
                cerr << "could not read -m, --mnp-ratio" << endl;
                exit(1);
            }
            break;
 
          case 'm':
            if (!convert(optarg, microsatellite_mutation_rate)) {
                cerr << "could not read -m, --microsat-rate" << endl;
                exit(1);
            }
            break;
 
          case 't':
            if (!convert(optarg, microsatellite_afs_alpha)) {
                cerr << "could not read -m, --microsatellite-afs-alpha" << endl;
                exit(1);
            }
            break;
 
          case 'j':
            if (!convert(optarg, microsatellite_len_alpha)) {
                cerr << "could not read -m, --microsatellite-len-alpha" << endl;
                exit(1);
            }
            break;
 
          case 'l':
            if (!convert(optarg, microsatellite_min_length)) {
                cerr << "could not read -l, --microsat-min-len" << endl;
                exit(1);
            }
            break;
 
          case 'p':
            if (!convert(optarg, ploidy)) {
                cerr << "could not read -p, --ploidy" << endl;
                exit(1);
            }
            break;

          case 'P':
            file_prefix = optarg;
            break;

          case 'S':
            sample_prefix = optarg;
            break;
 
          case 'n':
            if (!convert(optarg, population_size)) {
                cerr << "could not read -n, --population-size" << endl;
                exit(1);
            }
            sample_id_max_digits = strlen(optarg);
            break;

          case 'g':
            if (!convert(optarg, seed)) {
                cerr << "could not read -g, --random-seed" << endl;
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

    init_genrand(seed); // seed mt with current time

    string seqname;
    string sequence;  // holds sequence so we can process it

    FastaReference fr;
    fr.open(fastaFileName);

    string bases = "ATGC";

    vcf::VariantCallFile vcfFile;

    // write the VCF header
    stringstream headerss;
    headerss 
        << "##fileformat=VCFv4.1" << endl
        << "##fileDate=" << dateStr() << endl
        << "##source=mutatrix population genome simulator" << endl
        << "##seed=" << seed << endl
        << "##reference=" << fastaFileName << endl
        << "##phasing=true" << endl
        << "##commandline=" << command_line << endl
        << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate allele count\">" << endl
        << "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Type of each allele (snp, ins, del, mnp, complex)\">" << endl
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples at the site\">" << endl
        << "##INFO=<ID=NA,Number=1,Type=Integer,Description=\"Number of alternate alleles\">" << endl
        << "##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"Length of each alternate allele\">" << endl
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

    vector<string> samples;
    for (int i = 0; i < population_size; ++i) {
        stringstream sampless;
        sampless << sample_prefix << setfill('0') << setw(sample_id_max_digits) << i + 1; // one-based sample names
        samples.push_back(sampless.str());
        headerss << "\t" << sampless.str();
    }
    headerss << endl;

    // and set up our VCF output file
    string header = headerss.str();
    vcfFile.openForOutput(header);
    cout << vcfFile.header;

    int copies = ploidy * population_size;

    map<string, vector<SampleFastaFile*> > sequencesByRefseq;

    if (!dry_run) {
        for (FastaIndex::iterator s = fr.index->begin(); s != fr.index->end(); ++s) {

            FastaIndexEntry& indexEntry = s->second;
            seqname = indexEntry.name;

            vector<SampleFastaFile*>& sequences = sequencesByRefseq[seqname];
            for (int i = 0; i < population_size; ++i) {
                stringstream sname;
                sname << sample_prefix << setfill('0') << setw(sample_id_max_digits) << i + 1;
                string samplename = sname.str();
                for (int j = 0; j < ploidy; ++j) {
                    stringstream cname;
                    cname << j;
                    string chromname = cname.str();
                    string fullname = samplename + ":" + seqname + ":" + chromname;
                    string filename = file_prefix + fullname + ".fa";
                    //sequences.push_back(SampleFastaFile(filename, seqname));
                    sequences.push_back(new SampleFastaFile(filename, seqname));
                }
            }
        }
    }

    for (FastaIndex::iterator s = fr.index->begin(); s != fr.index->end(); ++s) {

        FastaIndexEntry& indexEntry = s->second;
        seqname = indexEntry.name;
        sequence = fr.getSequence(s->first);

        vector<SampleFastaFile*>& sequences = sequencesByRefseq[seqname];
        //sequences.resize(copies);
        
        long int pos = 0;
        long int microsatellite_end_pos = 0;
        while (pos < sequence.size() - 1) {

            string ref = sequence.substr(pos, 1); // by default, ref is just the current base

            // skip non-DNA sequence information
            if (!(ref == "A" || ref == "T" || ref == "C" || ref == "G")) {
                pos += ref.size();
                continue;
            }

            vector<Allele> alleles;

            // establish if we are in a repeat
            // and what motif is being repeated, how many times

            int len = 1;

            // get reference repeats
            // if we have a repeat, adjust the mutation rate
            // using length and direction-dependent
            // formula from "Likelihood-Based Estimation of Microsatellite Mutation Rates"
            // http://www.genetics.org/cgi/content/full/164/2/781#T1

            if (pos > microsatellite_end_pos) {

                map<string, int> repeats = repeatCounts(pos, sequence, repeat_size_max);
                if (repeats.size() > 0) {

                    string seq;
                    int repeat_count = 0;
                    // get the "biggest" repeat, the most likely ms allele at this site
                    for (map<string, int>::iterator r = repeats.begin(); r != repeats.end(); ++r) {
                        if (repeat_count < r->second) {
                            repeat_count = r->second;
                            seq = r->first;
                        }
                    }

                    int microsatellite_length = repeat_count * seq.size();

                    // record end of microsatellite so we don't generate more mutations until we pass it
                    microsatellite_end_pos = pos + microsatellite_length;

                    if (microsatellite_length > microsatellite_min_length
                            && genrand_real1() / copies 
                                < microsatellite_mutation_rate * repeat_count) {

                        // establish the relative rate of ins and del events
                        /*
                        long double repeatMutationDelProbability = microsatelliteDelProb(repeat_count);
                        long double repeatMutationInsProbability = microsatelliteInsProb(repeat_count);
                        long double indel_balance = 1;
                        if (repeatMutationInsProbability > repeatMutationDelProbability) {
                            indel_balance = repeatMutationInsProbability / repeatMutationDelProbability;
                        } else {
                            indel_balance = 1 - (repeatMutationInsProbability / repeatMutationDelProbability);
                        }
                        */
                        double indel_balance = 0.5;

                        // how many alleles at the site?

                        //int numalleles = min((int) floor(zetarandom(microsatellite_afs_alpha)), (int) ((double) repeat_count * indel_balance));
                        int numalleles = random_allele_frequency(repeat_count, microsatellite_afs_alpha);
                        //cout << "repeat_count: " << repeat_count << " numalleles: " << numalleles << endl;

                        map<int, bool> allele_lengths;
                        // lengths of the alleles
                        while (allele_lengths.size() < numalleles) {
                            int allele_length;
                            // TODO adjust length so that shorter events are more likely...
                            if (genrand_real1() > indel_balance) {
                                allele_length = -1 * min((int) floor(zetarandom(microsatellite_len_alpha)), repeat_count);
                            } else {
                                allele_length = min((int) floor(zetarandom(microsatellite_len_alpha)), repeat_count);
                            }
                            //cout << allele_length << endl;
                            map<int, bool>::iterator f = allele_lengths.find(allele_length);
                            if (f == allele_lengths.end()) {
                                allele_lengths[allele_length] = true;
                            }
                        }

                        // generate alleles
                        for (map<int, bool>::iterator f = allele_lengths.begin();
                                f != allele_lengths.end(); ++f) {

                            int allele_length = f->first;
                            int c = abs(f->first);
                            string alt = seq;

                            for (int i = 1; i < c; ++i)
                                alt += seq;

                            if (allele_length > 0) {
                                alleles.push_back(Allele(ref, ref + alt, "MICROSAT"));
                            } else {
                                alleles.push_back(Allele(ref + alt, ref, "MICROSAT"));
                            }
                        }
                        //cout << "alleles.size() == " << alleles.size() << endl;
                    }
                }
            }

            // snp case
            if (genrand_real1() / copies < snp_mutation_rate) {

                // make an alternate allele
                string alt = ref;
                while (alt == ref) {
                    alt = string(1, bases.at(genrand_int32() % 4));
                }

                if (genrand_real1() < mnp_ratio) {
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
                alleles.push_back(Allele(ref, alt));
            }

            // indel case
            if (genrand_real1() / copies < indel_mutation_rate) {
                // how many bp?
                len = (int) floor(zetarandom(indel_alpha));
                // guard against out-of-sequence indels
                if (pos + len < sequence.size() && len <= indel_max) {
                    if (genrand_int32() % 2 == 0) {
                    // deletion
                        alleles.push_back(Allele(sequence.substr(pos, 1 + len), sequence.substr(pos, 1)));
                    } else {
                        string alt;
                    // insertion?
                        // insert some random de novo bases
                        while (alt.length() < len + 1) {
                            alt += string(1, bases.at(genrand_int32() % 4));
                        }
                        alleles.push_back(Allele(ref, alt));
                    }
                } else {
                    // fall through
                }
            }

            // no mutation generated
            if (alleles.empty()) {
                for (int i = 0; i < copies; ++i) {
                    if (!dry_run) {
                        sequences.at(i)->write(ref);
                    }
                }
                pos += ref.size();
            } else {

                // TODO randomly distribute all the alleles throughout the population
                // generate allele frequencies for each
                // fun times...

                string genotype;

                vector<bool> alts;
                random_shuffle(alleles.begin(), alleles.end());

                vector<Allele> population_alleles;
                vector<Allele> present_alleles; // filtered for AFS > 0 in the sample
                
                // AFS simulation
                int remaining_copies = copies;
                while (remaining_copies > 0 && !alleles.empty()) {
                    Allele allele = alleles.back();
                    alleles.pop_back();
                    int allele_freq = random_allele_frequency(remaining_copies, afs_alpha);
                    if (allele_freq > 0) {
                        present_alleles.push_back(allele);
                        for (int i = 0; i < allele_freq; ++i) {
                            population_alleles.push_back(allele);
                        }
                        remaining_copies -= allele_freq;
                    }
                }

                // reference alleles take up the rest
                Allele reference_allele = Allele(ref, ref);
                for (int i = 0; i < remaining_copies; ++i) {
                    population_alleles.push_back(reference_allele);
                }
                //present_alleles.push_back(reference_allele);

                assert(population_alleles.size() == copies);

                // shuffle the alleles around the population
                random_shuffle(population_alleles.begin(), population_alleles.end());

                vcf::Variant var(vcfFile);
                var.sequenceName = seqname;
                var.position = pos + 1;
                var.quality = 99;
                var.id = ".";
                var.filter = ".";
                var.info["NS"].push_back(convert(population_size));
                var.info["NA"].push_back(convert(present_alleles.size()));
                var.format.push_back("GT");

                // establish the correct reference sequence and alternate allele set
                var.ref = ref;
                for (vector<Allele>::iterator a = present_alleles.begin(); a != present_alleles.end(); ++a) {
                    Allele& allele = *a;
                    if (allele.ref.size() > var.ref.size()) {
                        var.ref = allele.ref;
                    }
                }

                map<string, int> alleleIndexes;
                alleleIndexes[convert(reference_allele)] = 0; // XXX should we handle this differently, by adding the reference allele to present_alleles?
                int i = 1;
                for (vector<Allele>::iterator a = present_alleles.begin(); a != present_alleles.end(); ++a, ++i) {
                    Allele& allele = *a;
                    //cout << allele << " " << i << endl;
                    alleleIndexes[convert(allele)] = i;
                    //cout << allele << " " << i << endl;
                }

                //for (map<string, int>::iterator a = alleleIndexes.begin(); a != alleleIndexes.end(); ++a) {
                //    cout << a->first << " = " << a->second << endl;
                //}

                // now the reference allele is the largest possible, adjust the alt allele strings to reflect this
                // if we have indels, add the base before, set the position back one
                for (vector<Allele>::iterator a = present_alleles.begin(); a != present_alleles.end(); ++a) {
                    Allele& allele = *a;
                    string alleleStr = var.ref;
                    if (allele.ref.size() == allele.alt.size()) {
                        alleleStr.replace(0, allele.alt.size(), allele.alt);
                    } else {
                        alleleStr.replace(0, allele.ref.size(), allele.alt);
                    }
                    var.alt.push_back(alleleStr);
                }

                int j = 0;
                for (vector<string>::iterator s = samples.begin(); s != samples.end(); ++s, ++j) {
                    string& sample = *s;
                    vector<string> genotype;
                    // XXX hack, maybe this should get stored in another map for easier access?
                    for (int i = 0; i < ploidy; ++i) {
                        int l = (j * ploidy) + i;
                        //cout << l << " " << population_alleles.at(l) << " " << alleleIndexes[convert(population_alleles.at(l))] << endl;
                        genotype.push_back(convert(alleleIndexes[convert(population_alleles.at(l))]));
                    }
                    var.samples[sample]["GT"].push_back(join(genotype, "|"));
                    //cout << var.samples[sample]["GT"].front() << endl;
                }

                // write the genotype specs

                // now write out our sequence data (FASTA files), tabulate
                // allele frequency, and write some details to the VCF
                for (vector<Allele>::iterator a = present_alleles.begin(); a != present_alleles.end(); ++a) {

                    Allele& allele = *a;

                    vector<string> genotypes;
                    genotypes.resize(population_size);

                    int allele_freq = 0;

                    // obtain allele frequencies and output FASTA sequence data
                    // for each simulated sample
                    for (int j = 0; j < population_size; ++j) {
                        for (int i = 0; i < ploidy; ++i) {
                            int l = (j * ploidy) + i;
                            if (population_alleles.at(l) == allele) {
                                if (!dry_run) {
                                    sequences.at(l)->write(allele.alt);
                                }
                                ++allele_freq;
                            } else if (population_alleles.at(l) == reference_allele) {
                                if (!dry_run) {
                                    sequences.at(l)->write(allele.ref);
                                }
                            }
                        }
                    }

                    // set up the allele-specific INFO fields in the VCF record
                    var.info["AC"].push_back(convert(allele_freq));

                    int delta = allele.alt.size() - allele.ref.size();
                    if (delta == 0) {
                        if (allele.ref.size() == 1) {
                            var.info["TYPE"].push_back("snp");
                            var.info["LEN"].push_back(convert(allele.ref.size()));
                        } else {
                            var.info["TYPE"].push_back("mnp");;
                            var.info["LEN"].push_back(convert(allele.ref.size()));
                        }
                    } else if (delta > 0) {
                        var.info["TYPE"].push_back("ins");;
                        var.info["LEN"].push_back(convert(abs(delta)));
                    } else {
                        var.info["TYPE"].push_back("del");;
                        var.info["LEN"].push_back(convert(abs(delta)));
                    }
                    if (!allele.type.empty()) {
                        var.infoFlags[allele.type] = true;
                    }

                }

                // write the VCF record to stdout
                cout << var << endl;

                int largest_ref = 1; // enforce one pos
                for (vector<Allele>::iterator a = present_alleles.begin(); a != present_alleles.end(); ++a) {
                    if (a->ref.size() > largest_ref) {
                        largest_ref = a->ref.size();
                    }
                }

                pos += largest_ref; // step by the size of the last event
            }
        }
    }

    // close, clean up files
    for (map<string, vector<SampleFastaFile*> >::iterator s = sequencesByRefseq.begin(); s != sequencesByRefseq.end(); ++s) {
        vector<SampleFastaFile*>& files = s->second;
        for (vector<SampleFastaFile*>::iterator f = files.begin(); f != files.end(); ++f) {
            delete *f;
        }
        files.clear();
    }

    return 0;

}

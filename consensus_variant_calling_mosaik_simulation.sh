#!/bin/bash


if [ $# -lt 5 ]; then
    echo usage: $0 [ref] [num samples] [coverage] [readlength] [fragment length]
    exit
fi

reference=$1
n=$2 # number of samples
coverage=$3
readlength=$4
mfl=$5
mosaik_reference=`echo $reference | sed -e "s/.fa/.mbr/"`
answers=answers.vcf.gz


echo
echo generating test data
echo

# random seed for mutatrix
seed=$(od -An -N2 -i /dev/random)
mutatrix -g $seed -S sample -p 2 -n $n $reference | bgziptabix $answers
samples=`zcat $answers | grep ^#CHROM | cut -f $n- | sed -e "s/\t/\n/g" `
samples=$(for i in $(seq -w $n); do echo sample$i; done)
echo $samples

# calculate coverage
nbp=$(cat $reference | grep -v ">" | sed "s/\n//g" | wc -c)
nreads=$(($nbp * $coverage / $readlength))

# generate reads given no haplotype mutations

echo
echo simulating and aligning...
echo

MosaikBuild -fr $reference -oa $mosaik_reference
refname=`head -1 $reference | sed -e "s/>//"`
for sample in $samples; do cat $sample:$refname:*.fa >$sample.fa; done
for sample in $samples;
do
    seed=$(od -An -N2 -i /dev/random)
    mason illumina -n $readlength -N $nreads -rnp $sample -s $seed -sq -mp -ll $mfl -hs 0 -hi 0 $sample.fa
done
for sample in $samples; do MosaikBuild -st illumina -sam $sample -id $sample -q $sample.fa_1.fastq -q2 $sample.fa_2.fastq -out $sample.mra; done
for sample in $samples; do echo $sample; MosaikAligner -p 8 -mm 16 -in $sample.mra -out $sample -ia $mosaik_reference; done
for sample in $samples
do
    bamtools sort -in $sample.bam -out $sample.temp.bam
    bamtools resolve -twoPass -in $sample.temp.bam -out $sample.sorted.bam 
    rm $sample.temp.bam
done
for file in *sorted.bam; do bamtools index -in $file; done

echo
echo calling...
echo

results=results.vcf
answers_primitives=answers.primitives.vcf.gz

vcfallelicprimitives $answers | bgziptabix $answers_primitives


sets=""

files=""
for bam in *sorted.bam; do files=$files" -in $bam "; done

echo
echo freebayes
time freebayes -f $reference *sorted.bam >freebayes.$results 2>freebayes.timing &
sets="$sets freebayes.$results"

echo
echo freebayes k20
time freebayes -f $reference --max-complex-gap 20 *sorted.bam >freebayes.k20.$results 2>freebayes.timing &
sets="$sets freebayes.k20.$results"

# samtools is too slow to test on >100 samples, but here goes
echo
echo samtools
time samtoolsbam2vcf $reference *sorted.bam >samtools.$results 2>samtools.timing &
sets="$sets samtools.$results"

echo
echo gatk
time gatkbam2vcf $reference *sorted.bam | grep -v "^INFO" | grep -v "^WARN" >gatk.$results 2>gatk.timing &
sets="$sets gatk.$results"

echo
echo gatk open
files=""
for bam in *sorted.bam; do files=$files" -I $bam "; done
time java -jar /share/home/erik/src/GenomeAnalysisTKLite-2.2-15-g4828906/GenomeAnalysisTKLite.jar \
    -T UnifiedGenotyper -stand_emit_conf 1 -glm BOTH -R $reference $files | grep -v "^INFO" | grep -v "^WARN" >gatk.open.$results 2>gatk.open.timing &
sets="$sets gatk.open.$results"

echo
echo sets are $sets

echo
echo waiting
wait  # awaits completion of variant detection


echo -e "set\tthreshold\tnum_snps\tfalse_positive_snps\tfalse_negative_snps\tnum_indels\tfalse_positive_indels\tfalse_negative_indels" >results.roc.tsv

for results in $sets;
do

    set=$(basename $results .results.vcf)

    echo
    echo comparing $results from $set
    echo

    vcfallelicprimitives $results | vcfroc -r $reference -t $answers_primitives | grep -v threshold | sed "s/^/$set\t/" >>results.roc.tsv

done

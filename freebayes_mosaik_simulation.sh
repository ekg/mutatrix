#!/bin/bash


n=10
reference=chr20_bit.fa 
mosaik_reference=`echo $reference | sed -e "s/.fa/.mbr/"`
answers=answers.vcf

readlength=70
nreads=20000
mfl=500

echo
echo generating test data
echo
#mutatrix -r 0.0001 -a 2 -S sample -p 2 -n 10 $reference | tee $answers
mutatrix -S sample -p 2 -n 10 $reference | tee $answers
samples=`cat $answers | grep ^#CHROM | cut -f 10- | sed -e "s/\t/\n/g" `

echo
echo aligning...
echo
MosaikBuild -fr $reference -oa $mosaik_reference
refname=`head -1 $reference | sed -e "s/>//"`
for sample in $samples; do cat $sample:$refname:*.fa >$sample.fa; done
for sample in $samples; do wgsim -1 $readlength -2 $readlength -R 0.02 -r 0 -d $mfl -N $nreads $sample.fa  $sample.read1.fa $sample.read2.fa; done
#for sample in $samples; do art_illumina_q2e $sample.fa $sample.reads 70 10; done
for sample in $samples; do MosaikBuild -st illumina -sam $sample -id $sample -q $sample.read1.fa -q2 $sample.read2.fa -out $sample.mra; done
for sample in $samples; do echo $sample; MosaikAligner -p 8 -mm 16 -in $sample.mra -out $sample -ia $mosaik_reference; done
for sample in $samples; do bamtools sort -in $sample.bam -out $sample.sorted.bam; done
for file in *sorted.bam; do bamtools index -in $file; done

echo
echo calling...
echo

results=results.vcf
#bamtools filter $(for file in *sorted.bam; do echo $file; done | sed -e "s/^/ -in /g") \
#    | bamleftalign -f $reference \
#    | samtools fillmd -Aru - $reference 2>/dev/null \
#    | freebayes --no-filters -d -f $reference --stdin \
#    | tee $results.unsorted
freebayes -C 2 --no-filters --left-align-indels -f $reference *sorted.bam >$results

# sort the vcf file

echo
echo comparing...
echo
#bgzip -f $answers
#tabix -fp vcf $answers.gz
#bgzip -f $results
#tabix -fp vcf $results.gz

#echo Q >= 0
#vcf-compare $answers.gz $results.gz

for Q in 0 1 10 20 30 40
do
    echo
    echo "Q >= $Q"
    echo
    echo "records only in answers.vcf:"
    cat results.vcf | vcffilter -f "QUAL > $Q" | vcfCTools intersect -q a -i answers.vcf -i - | grep -v ^# | wc -l
    echo "records only in results.vcf:"
    cat results.vcf | vcffilter -f "QUAL > $Q" | vcfCTools intersect -q a -i - -i answers.vcf | grep -v ^# | wc -l
    echo
done

#!/bin/bash


n=10
reference=chr20_bit.fa 
mosaik_reference=`echo $reference | sed -e "s/.fa/.mbr/"`
answers=answers.vcf

echo
echo generating test data
echo
mutatrix -r 0.0001 -a 2 -S sample -p 2 -n 10 $reference | tee $answers
samples=`cat $answers | grep ^#CHROM | cut -f 10- | sed -e "s/\t/\n/g" `

echo
echo aligning...
echo
MosaikBuild -fr $reference -oa $mosaik_reference
refname=`head -1 $reference | sed -e "s/>//"`
for sample in $samples; do cat $sample:$refname:*.fa >$sample.fa; done
for sample in $samples; do wgsim -1 70 -2 70 -R 0.02 -r 0 -N 20000 $sample.fa  $sample.read1.fa $sample.read2.fa; done
#for sample in $samples; do art_illumina_q2e $sample.fa $sample.reads 70 10; done
for sample in $samples; do MosaikBuild -st illumina -sam $sample -id $sample -q $sample.read1.fa -q2 $sample.read2.fa -out $sample.mra; done
for sample in $samples; do echo $sample; MosaikAligner -p 8 -mm 16 -in $sample.mra -out $sample -ia $mosaik_reference; done
for sample in $samples; do bamtools sort -in $sample.bam -out $sample.sorted.bam; done
for file in *sorted.bam; do bamtools index -in $file; done

echo
echo calling...
echo

results=results.vcf
bamtools filter $(for file in *sorted.bam; do echo $file; done | sed -e "s/^/ -in /g") \
    | bamleftalign -f $reference \
    | samtools fillmd -Aru - $reference 2>/dev/null \
    | freebayes -U 4 -e 1 -C 2 --no-filters --mnps --indels -f $reference --stdin \
    | tee $results.unsorted

# sort the vcf file
vcf_header $results.unsorted >$results
cat $results.unsorted | grep -v ^# | sort -k1,1d -k2,2n >>$results

echo
echo comparing...
echo
bgzip -f $answers
tabix -fp vcf $answers.gz
bgzip -c $results >$results.gz
tabix -fp vcf $results.gz

echo Q >= 0
compare-vcf $answers.gz $results.gz

for Q in 10 20 30 40
do
    echo
    echo "Q >= $Q"
    echo
    resultsfilt=$results.$Q.vcf
    cat $results | vcfqualfilter -c $Q >$resultsfilt && bgzip -c $resultsfilt >$resultsfilt.gz && tabix -fp vcf $resultsfilt.gz
    compare-vcf -w 4 $answers.gz $resultsfilt.gz
    echo
done

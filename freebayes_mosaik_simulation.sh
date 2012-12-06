#!/bin/bash


n=1
reference=chr20_bit.fa 
mosaik_reference=`echo $reference | sed -e "s/.fa/.mbr/"`
answers=answers.vcf.gz

readlength=70
nreads=20000
mfl=500

echo
echo generating test data
echo
#mutatrix -r 0.0001 -a 2 -S sample -p 2 -n 10 $reference | tee $answers
mutatrix -S sample -p 2 -n 10 $reference | bgziptabix $answers
samples=`zcat $answers | grep ^#CHROM | cut -f 10- | sed -e "s/\t/\n/g" `

echo
echo aligning...
echo
MosaikBuild -fr $reference -oa $mosaik_reference
refname=`head -1 $reference | sed -e "s/>//"`
for sample in $samples; do cat $sample:$refname:*.fa >$sample.fa; done
for sample in $samples; do wgsim -1 $readlength -2 $readlength -e 0.02 -R 0.01 -r 0 -d $mfl -N $nreads $sample.fa  $sample.read1.fa $sample.read2.fa; done
#for sample in $samples; do wgsim -1 $readlength -2 $readlength -e 0.05 -R 0.05 -r 0 -d $mfl -N $nreads $sample.fa  $sample.read1.fa $sample.read2.fa; done
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
answers_primitives=answers.primitives.vcf.gz
vcfallelicprimitives $answers | bgziptabix $answers_primitives
#freebayes.unstable -C 2 --max-complex-gap 40 --no-filters --left-align-indels -f $reference --haplotype-basis-alleles $answers_primitives *sorted.bam >$results
freebayes.unstable -C 2 --left-align-indels -f $reference *sorted.bam >$results

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

#vcf2bed.py <answers.vcf >answers.bed

echo -e QUAL\\tnum_sites\\tfalse_positive_sites\\tfalse_negative_sites\\tnum_snps\\tfalse_positive_snps\\tfalse_negative_snps\\tnum_mnps\\tfalse_positive_mnps\\tfalse_negative_mnps\\tnum_indels\\tfalse_positive_indels\\tfalse_negative_indels >results.roc.tsv
cat results.roc.tsv

last=results.vcf
levels=$(seq 0 1 100)

for Q in $levels 
do
    cat $last | vcffilter -f "QUAL > $Q" | vcfallelicprimitives >results.$Q.vcf
    last=results.$Q.vcf

    vcfintersect -r $reference -v -i results.$Q.vcf $answers_primitives | vcfstats >false_negatives.$Q.stats
    vcfintersect -r $reference -v -i $answers_primitives results.$Q.vcf | vcfstats >false_positives.$Q.stats
    vcfstats results.$Q.vcf >calls.$Q.stats

    fn=$(cat false_negatives.$Q.stats | grep "total variant sites:" | cut -f 2)
    fp=$(cat false_positives.$Q.stats | grep "total variant sites:" | cut -f 2)
    nc=$(cat calls.$Q.stats | grep "total variant sites:" | cut -f 2)

    fns=$(cat false_negatives.$Q.stats | grep "^snps:" | cut -f 2)
    fps=$(cat false_positives.$Q.stats | grep "^snps:" | cut -f 2)
    ncs=$(cat calls.$Q.stats | grep "^snps:" | cut -f 2)

    fnm=$(cat false_negatives.$Q.stats | grep "^mnps:" | cut -f 2)
    fpm=$(cat false_positives.$Q.stats | grep "^mnps:" | cut -f 2)
    ncm=$(cat calls.$Q.stats | grep "^mnps:" | cut -f 2)

    fni=$(cat false_negatives.$Q.stats | grep "^indels:" | cut -f 2)
    fpi=$(cat false_positives.$Q.stats | grep "^indels:" | cut -f 2)
    nci=$(cat calls.$Q.stats | grep "^indels:" | cut -f 2)

    rm false_negatives.$Q.stats false_positives.$Q.stats calls.$Q.stats


    echo -e $Q\\t$nc\\t$fp\\t$fn\\t$ncs\\t$fps\\t$fns\\t$ncm\\t$fpm\\t$fnm\\t$nci\\t$fpi\\t$fni >>results.roc.tsv
    echo -e $Q\\t$nc\\t$fp\\t$fn\\t$ncs\\t$fps\\t$fns\\t$ncm\\t$fpm\\t$fnm\\t$nci\\t$fpi\\t$fni
done

for Q in $levels
do
    rm results.$Q.vcf
done

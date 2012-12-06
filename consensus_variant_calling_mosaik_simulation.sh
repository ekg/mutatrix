#!/bin/bash



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
#bamtools filter $(for file in *sorted.bam; do echo $file; done | sed -e "s/^/ -in /g") \
#    | bamleftalign -f $reference \
#    | samtools fillmd -Aru - $reference 2>/dev/null \
#    | freebayes --no-filters -d -f $reference --stdin \
#    | tee $results.unsorted
answers_primitives=answers.primitives.vcf.gz
vcfallelicprimitives $answers | bgziptabix $answers_primitives
#freebayes.unstable -C 2 --max-complex-gap 40 --no-filters --left-align-indels -f $reference --haplotype-basis-alleles $answers_primitives *sorted.bam >$results

files=""
for bam in *sorted.bam; do files=$files" -in $bam "; done

echo freebayes
time ( bamtools merge $files \
    | bamleftalign -f $reference \
    | samtools calmd -EAru - $reference 2>/dev/null \
    | freebayes -C 2 --stdin -f $reference >freebayes.default.$results ) &

echo freebayes maf0.3
time ( bamtools merge $files \
    | bamleftalign -f $reference \
    | samtools calmd -EAru - $reference 2>/dev/null \
    | freebayes -C 2 -F 0.3 --stdin -f $reference >freebayes.default.maf0.3.$results ) &

echo freebayes ogap
time ( bamtools merge $files \
    | ogap -z -R 25 -C 20 -Q 20 -S 0 -f $reference \
    | bamleftalign -f $reference \
    | samtools calmd -EAru - $reference 2>/dev/null \
    | freebayes -C 2 --stdin -f $reference >freebayes.default.ogap.$results ) &

echo freebayes local assembly
time ( bamtools merge $files \
    | ogap -z -R 25 -C 20 -Q 20 -S 0 -f $reference \
    | bamleftalign -f $reference \
    | samtools calmd -EAru - $reference 2>/dev/null \
    | freebayes -C 2 --stdin --max-complex-gap 30 -f $reference >freebayes.local_assembly.$results ) &

echo
echo samtools
time samtoolsbam2vcf $reference *sorted.bam >samtools.$results &

echo
echo gatk
time gatkbam2vcf $reference *sorted.bam | grep -v "^INFO" | grep -v "^WARN" >gatk.$results &

echo
echo waiting
wait  # awaits completion of variant detection


# sort the vcf file

echo -e set\\tQUAL\\tnum_sites\\tfalse_positive_sites\\tfalse_negative_sites\\tnum_snps\\tfalse_positive_snps\\tfalse_negative_snps\\tnum_mnps\\tfalse_positive_mnps\\tfalse_negative_mnps\\tnum_indels\\tfalse_positive_indels\\tfalse_negative_indels >results.roc.tsv
cat results.roc.tsv


for results in freebayes.default.$results freebayes.default.maf0.3.$results freebayes.default.ogap.$results freebayes.local_assembly.$results samtools.$results gatk.$results;
do

    set=$(basename $results .results.vcf)

    echo
    echo comparing $results from $set
    echo

    last=$results
    levels=$(seq 0 1 100)

    for Q in $levels 
    do
	cat $last | vcffilter -f "QUAL > $Q" | vcfallelicprimitives >$results.$Q.vcf
	last=$results.$Q.vcf

	vcfintersect -r $reference -v -i $results.$Q.vcf $answers_primitives | vcfstats >false_negatives.$Q.stats
	vcfintersect -r $reference -v -i $answers_primitives $results.$Q.vcf | vcfstats >false_positives.$Q.stats
	vcfstats $results.$Q.vcf >calls.$Q.stats

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


	echo -e $set\\t$Q\\t$nc\\t$fp\\t$fn\\t$ncs\\t$fps\\t$fns\\t$ncm\\t$fpm\\t$fnm\\t$nci\\t$fpi\\t$fni >>results.roc.tsv
	echo -e $set\\t$Q\\t$nc\\t$fp\\t$fn\\t$ncs\\t$fps\\t$fns\\t$ncm\\t$fpm\\t$fnm\\t$nci\\t$fpi\\t$fni
    done

    for Q in $levels
    do
	rm $results.$Q.vcf
    done

done

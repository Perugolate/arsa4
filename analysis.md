# Table of contents

- [Dependencies](#dependencies)
- [Get data, do QC](#get-data-do-qc)
- [Prepare reference](#prepare-reference)
- [Call SNPs](#call-snps)


## Dependencies

* prokka >=1.12-beta (and dependencies)
* snippy >=2.9
* Perl >= 5.6
* BioPerl >= 1.6
* bwa mem >= 0.7.12 
* samtools >= 1.1
* GNU parallel > 2013xxxx
* freebayes >= 1.0.1 
* freebayes sripts (freebayes-parallel, fasta_generate_regions.py)
* vcflib (vcffilter, vcfstreamsort, vcfuniq, vcffirstheader)
* vcftools (vcf-consensus)
* snpEff >= 4.1

```bash
source paths
```

## Get data, do QC

```bash
cd $OMDAT
basemount $BASESPACE
cp $BASESPACE/Projects/2015_29_JR_bacterial_genomes/Samples/*/Files/Data/Intensities/BaseCalls/*fastq.gz .
fastqc *fastq.gz
mkdir qc && mv *zip qc/ && mv *html qc/
```

## Prepare reference

Re-use 8325-4 reference genome from ARSA3. Three prophages have been excised after Baek et al. 2013.
```bash
cd $ARSA4
mkdir ref && cd ref
cp $CAROZ/jlal513/flat/phage1/jla513_unpol/jla513_unpol.fna .
prokka --cpus 12 jla513_unpol.fna --prefix jla513_unpol
```

Call SNPs in ancestor and correct reference.
```bash
cd $ARSA4/ref
snippy --cpus 12 --outdir 76_S74 --ref jla513_unpol/jla513_unpol.gbk --rgid 76_S74 --prefix 76_S74 --pe1 $OMDAT/76_S74_L001_R1_001.fastq.gz --pe2 $OMDAT/76_S74_L001_R2_001.fastq.gz &> 76_S74.log
bcftools consensus -f jla513_unpol.fna  76_S74/76_S74.vcf.gz  -o sh1000_pol_1.fna
prokka --cpus 12 sh1000_pol_1.fna --prefix sh1000_pol_1
```

Call SNPs in present in unselected controls.
```bash
grep -e WT -e G $ARSA4/metadata.csv | cut -f2 -d"," | sed -e 's/$/_S/g' -e 's/^/\^/g' | grep -f - $ARSA4/samples.txt > controls.txt
for i in `cat controls.txt`
do
  snippy --cpus 12 --outdir $i --ref sh1000_pol_1.gbk --rgid $i --prefix $i --pe1 $OMDAT/${i}_L001_R1_001.fastq.gz --pe2 $OMDAT/${i}_L001_R2_001.fastq.gz &> $i.log
done
```

Determine core SNPs present in all unselected controls.
```bash
snippy-core --prefix control_core *_S*/
```
There are many common SNPs in the controls but none of them are present in every single control strain. 

## Call SNPs

Call SNPs in selected lines.
```bash
cd $OMDAT
ls *gz | cut -f1-2 -d"_" | sort -u > ../an_b1/list
cd $ARSA4
for i in `cat list`
do
  snippy --cpus 12 --outdir $i --ref sh1000_pol_1.gbk --rgid $i --prefix $i  --pe1 $OMDAT/${i}_L001_R1_001.fastq.gz --pe2 $OMDAT/${i}_L001_R2_001.fastq.gz &> $i.log
done
snippy-core --prefix coreAll *_S*/
```


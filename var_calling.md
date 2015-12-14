# Table of contents

- [Dependencies](#dependencies)
- [Get data, do QC](#get-data-do-qc)
- [Prepare reference](#prepare-reference)
- [Call SNPs](#call-snps)
 * [Call SNPs with `freebayes` using `snippy`
](#call-snps-with-freebayes-using-snippy)
 * [Call variants with `breseq`](#call-variants-with-breseq)
- [CNV analysis](#cnv-analysis)

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

Re-use 8325 reference genome from ARSA3. Three prophages have been excised after Baek et al. 2013.
```bash
cd $ARSA4
mkdir ref && cd ref
cp $CAROZ/jlal513/flat/phage1/no_phi.fna .
prokka --cpus 12 no_phi.fna --prefix no_phi
```

Call SNPs in ancestor and correct reference.
```bash
cd $ARSA4/ref
snippy --cpus 12 --outdir 76_S74 --ref no_phi/no_phi.gbk --rgid 76_S74 --prefix 76_S74 --pe1 $OMDAT/76_S74_L001_R1_001.fastq.gz --pe2 $OMDAT/76_S74_L001_R2_001.fastq.gz &> 76_S74.log
bcftools consensus -f no_phi.fna  76_S74/76_S74.vcf.gz  -o sh1000_pol_1.fna
prokka --cpus 12 sh1000_pol_1.fna --prefix sh1000_pol_1
```

The resultsing polished reference still had some remaining variants. Mostly medium deletions that were not called by `snippy`. These were iteratively identified, manually corrected, and re-checked using `breseq`. The resulting polished reference is `sh1000_pol_6.fna`.

Call SNPs in present in unselected controls.
```bash
grep -e WT -e G $ARSA4/metadata.csv | cut -f2 -d"," | sed -e 's/$/_S/g' -e 's/^/\^/g' | grep -f - $ARSA4/samples.txt > controls.txt
for i in `cat controls.txt`
do
  snippy --cpus 12 --outdir $i --ref sh1000_pol_6.gbk --rgid $i --prefix $i --pe1 $OMDAT/${i}_L001_R1_001.fastq.gz --pe2 $OMDAT/${i}_L001_R2_001.fastq.gz &> $i.log
done
```

Determine core SNPs present in all unselected controls.
```bash
snippy-core --prefix control_core *_S*/
```
This showed that there are no SNPs that are present in every single control strain. There are SNPs which segregate by the biological replicate line.

## Call SNPs

### Call SNPs with `freebayes` using `snippy`
```bash
cd $ARSA4
grep T1 $ARSA4/metadata.csv | cut -f2 -d"," | sed -e 's/$/_S/g' -e 's/^/\^/g' | grep -f - $ARSA4/samples.txt > selected.txt
for i in `cat selected.txt`
do
  snippy --cpus 12 --outdir $i --ref sh1000_pol_6.gbk --rgid $i --prefix $i  --pe1 $OMDAT/${i}_L001_R1_001.fastq.gz --pe2 $OMDAT/${i}_L001_R2_001.fastq.gz &> $i.log
done
snippy-core --prefix coreAll *_S*/
```

### Call variants with `breseq`

`breseq` is slow so use HPC. Template job SLURM file:

```bash
#! /bin/bash
#SBATCH --job-name=template_breseq
#SBATCH --mem=1024
#SBATCH --cpus-per-task=4
#SBATCH --time=00:20:00
cd /scratch/$USER/omdat/
breseq -r sh1000_pol_6.gbk template_L001_R1_001.fastq template_L001_R2_001.fastq -j 4 -n template -o template
```

Note 20 minutes is cutting it fine but there were no failures.

Generate job for each sample.
```bash
for i in `cat samples.txt`
do
  sed "s/template/$i/g" template.sh > $i.sh
done
```

```bash
cd $ARSA4
mkdir breseq && breseq
for i in `cat ../samples.txt`
do
  scp $SOROBAN:/scratch/$USER/olga90/$i/output/index.html $i.index.html
done
```

## CNV analysis

Prepare alignments.
```bash
cd $ARSA4
mkdir cnv && cd cnv
cp ../ref/sh1000_pol_6/sh1000_pol_6.fna .
bwa index sh1000_pol_6.fna
```

Align reads to reference, sort, and index.
```bash
for i in `cat ../samples.txt`
do
  bwa mem -t 12 -R "@RG\tID:${i}\tSM:${i}" sh1000_pol_6.fna $OMDAT/${i}_L001_R1_001.fastq.gz $OMDAT/${i}_L001_R2_001.fastq.gz > $i.sam
  picard-tools SortSam INPUT=$i.sam OUTPUT=$i.bam SORT_ORDER=coordinate
  picard-tools BuildBamIndex INPUT=$i.bam
done
```

Analyse CNVs in `R`.
```R
library("cn.mops")
library("magrittr")
BAMFiles <- list.files(pattern=".bam$")
bamDataRanges <- getReadCountsFromBAM(BAMFiles, mode="paired")
res <- haplocn.mops(bamDataRanges)
res <- calcIntegerCopyNumbers(res)
plot(res, which=1)
```


# Table of contents

- [Dependencies](#dependencies)
- [Get data, do QC](#get-data-do-qc)
- [Prepare reference](#prepare-reference)
- [Call SNPs](#call-snps)
 * [Call SNPs with `freebayes` using `snippy`](#call-snps-with-freebayes-using-snippy)
 * [Call variants with `breseq`](#call-variants-with-breseq)
- [CNV analysis](#cnv-analysis)

## Dependencies

* prokka >=1.12-beta (and dependencies)
* snippy >=3.0
* Perl >= 5.6
* BioPerl >= 1.6
* bwa mem >= 0.7.12 
* samtools >= 1.3
* GNU parallel > 2013xxxx
* freebayes >= 1.0.1 
* freebayes scripts (freebayes-parallel, fasta_generate_regions.py)
* vcflib (vcffilter, vcfstreamsort, vcfuniq, vcffirstheader)
* vcftools (vcf-consensus)
* snpEff >= 4.1

```sh
source paths
```

## Get data, do QC

**Sort this section out - merging the reads by sample**

```sh
cd $OMDAT
basemount $BASESPACE
cp $BASESPACE/Projects/2015_29_JR_bacterial_genomes/Samples/*/Files/Data/Intensities/BaseCalls/*fastq.gz .
fastqc *fastq.gz
cp $BASESPACE/Projects/2015_29_JR_bacterial_genomes_II/Samples/*/Files/Data/Intensities/BaseCalls/*fastq.gz .
mkdir qc && mv *zip qc/ && mv *html qc/
```

## Prepare reference

Re-use 8325 reference genome from ARSA3. Three prophages have been excised after Baek et al. 2013.
```sh
cd $ARSA4
mkdir ref && cd ref
cp $CAROZ/jlal513/flat/phage1/no_phi.fna .
prokka --cpus 12 no_phi.fna --prefix no_phi
```

Call SNPs in ancestor and correct reference.
```sh
cd $ARSA4/ref
snippy --cpus 12 --outdir 76_S74 --ref no_phi/no_phi.gbk --rgid 76_S74 --prefix 76_S74 --pe1 $OMDAT/76_S74_L001_R1_001.fastq.gz --pe2 $OMDAT/76_S74_L001_R2_001.fastq.gz &> 76_S74.log
bcftools consensus -f no_phi.fna  76_S74/76_S74.vcf.gz  -o sh1000_pol_1.fna
prokka --cpus 12 sh1000_pol_1.fna --prefix sh1000_pol_1
```

The resultsing polished reference still had some remaining variants. Mostly medium deletions that were not called by `snippy`. These were iteratively identified, manually corrected, and re-checked using `breseq`. The resulting polished reference is `sh1000_pol_6.fna`.

Call SNPs in present in unselected controls.
```sh
grep -e WT -e G $ARSA4/metadata.csv | cut -f2 -d"," | sed -e 's/$/_S/g' -e 's/^/\^/g' | grep -f - $ARSA4/samples.txt > controls.txt
for i in `cat controls.txt`
do
  snippy --cpus 12 --outdir $i --ref sh1000_pol_6.gbk --rgid $i --prefix $i --pe1 $OMDAT/${i}_L001_R1_001.fastq.gz --pe2 $OMDAT/${i}_L001_R2_001.fastq.gz &> $i.log
done
```

Determine core SNPs present in all unselected controls.
```sh
snippy-core --prefix control_core *_S*/
```
This showed that there are no SNPs that are present in every single control strain. There are SNPs which segregate by the biological replicate line.

## Trimming

**No need to do trimming/calling of controls and selected separately - this is just a hangover from previous version**

```sh
for i in `cat $ARSA4/controls.txt`
do
  java -jar $TRIMO PE -baseout $OMDAT/trimmed/${i}.fq.gz -threads 12 -phred33 $OMDAT/${i}_$R1 $OMDAT/${i}_$R2 ILLUMINACLIP:/home/paul/src/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
done
```

```sh
for i in `cat $ARSA4/selected.txt`
do
  java -jar $TRIMO PE -baseout $OMDAT/trimmed/${i}.fq.gz -threads 12 -phred33 $OMDAT/${i}_$R1 $OMDAT/${i}_$R2 ILLUMINACLIP:/home/paul/src/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
done
```

Note `snippy` does not allow use of both paired and unpaired reads.

## Call SNPs

### Call SNPs with `freebayes` using `snippy`
```sh
for i in `cat $ARSA4/controls.txt`
do
  snippy --cpus 12 --outdir $OMDAT/snippy/$i --ref $ARSA4/sh1000_pol_6.gbk --rgid $i --prefix $i --pe1 $OMDAT/trimmed/${i}_1P.fq.gz --pe2 $OMDAT/trimmed/${i}_2P.fq.gz &> $OMDAT/snippy/$i.log
done
# selected strains
cd $ARSA4
grep -e "T1" $ARSA4/metadata.csv | cut -f2 -d"," | sed -e 's/$/_S/g' -e 's/^/\^/g' | grep -f - $ARSA4/samples.txt > selected.txt
for i in `cat $ARSA4/selected.txt`
do
  snippy --cpus 12 --outdir $OMDAT/snippy/$i --ref $ARSA4/sh1000_pol_6.gbk --rgid $i --prefix $i --pe1 $OMDAT/trimmed/${i}_1P.fq.gz --pe2 $OMDAT/trimmed/${i}_2P.fq.gz &> $OMDAT/snippy/$i.log
done
snippy-core --prefix coreAll *_S*/
snippy-core --prefix coreAll *_S*/
```

Note, coverage is too low in 63_S61 (<59% at 10X).

### Call variants with `breseq`

`breseq` is slow so use HPC. Template job SLURM file:

```sh
#! /bin/sh
#SBATCH --job-name=template_breseq
#SBATCH --mem=1024
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
cd /scratch/$USER/omdat/
breseq -r sh1000_pol_6.gbk template_1P.fq.gz template_2P.fq.gz -j 4 -n template -o template
```

Note 20 minutes is cutting it fine but there were no failures.

Generate job for each sample.
```sh
for i in `cat samples.txt`
do
  sed "s/template/$i/g" template.sh > $i.sh
done
```

```sh
cd $ARSA4
mkdir breseq && breseq
for i in `cat ../samples.txt`
do
  scp $SOROBAN:/scratch/$USER/olga90/$i/output/index.html $i.index.html
done
```

## CNV analysis

Prepare alignments.
```sh
cd $ARSA4
mkdir cnv && cd cnv
cp $ARSA4/sh1000_pol_6.fa .
bwa index sh1000_pol_6.fa
```

Align reads to reference, sort, and index.
```sh
for i in `cat $ARSA4/samples.txt`
do
  bwa mem -t 12 -R "@RG\tID:${i}\tSM:${i}" sh1000_pol_6.fa $OMDAT/trimmed/${i}_1P.fq.gz $OMDAT/trimmed/${i}_2P.fq.gz > $i.sam
  picard-tools SortSam INPUT=$i.sam OUTPUT=$i.bam SORT_ORDER=coordinate
  picard-tools BuildBamIndex INPUT=$i.bam
done
```

Analyse CNVs in `R`.

This seems to miss some glaring deletions (unless they are mapping artefacts).
```R
library("cn.mops")
library("magrittr")
BAMFiles <- list.files(pattern = ".bam$")
bamDataRanges <- getReadCountsFromBAM(BAMFiles, mode = "paired", WL = 1000)
res <- haplocn.mops(bamDataRanges)
res <- calcIntegerCopyNumbers(res)
cnvs(res)
plot(res, which = 1)
```


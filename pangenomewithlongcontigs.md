### Pangenome with only contigs > 1M
I started with two strains.

*flaviav@penguin2:/home/flaviav/data/BXD/BXD005+BXD006_1M*

A. I extracted only the contigs greater than 1 million
```shell
awk -v n=1000000 '/^>/{ if(l>n) print b; b=$0;l=0;next }{l+=length;b=b ORS $0}END{if(l>n) print b }'  BXD005_supernova_changeid.fasta > BXD005_len1M.fa
awk -v n=1000000 '/^>/{ if(l>n) print b; b=$0;l=0;next }{l+=length;b=b ORS $0}END{if(l>n) print b }'  BXD006_supernova_changeid.fasta > BXD006_len1M.fa
```

#### 1) I mapped two strains against reference genome
 
*flaviav@penguin2:/home/flaviav/data/BXD/BXD005+BXD006_1M*

```shell
ls assemblies | grep '_len1M.fa$' | cut -f 1-2 -d . | sort -V | uniq >haps2strains.list
ref=assemblies/UCSC_mm10_changeid.fa
for hap in $(cat haps2strains.list);
do
    in=assemblies/$hap
    out=alignments/$hap.vs.ref.paf
    wfmash -t 16 -m -N -p 90 $ref $in >$out
done
```
#### 2) I split the contigs by chromosomes
```shell
(seq 19; echo X; echo Y) | while read i; do awk '$6 ~ "chr'$i'$"' $(ls alignments/*.vs.ref.paf | sort -V) | cut -f 1 | sort -V >parts/chr$i.contigs; done
```
####  3) I extracted only contigs of chr19

*flaviav@penguin2:/home/flaviav/data/BXD/BXD005+BXD006_1M/whole*

```shell

cat BXD005_len1M.fa BXD006_len1M.fa > BXD005+BXD006_len1M.fa
samtools faidx  BXD001+BXD002_len500kb.fa $(cat /home/flaviav/data/BXD/BXD005+BXD006_1M/parts/chr19.contigs) >chr19_BXD005+BXD006len1M.pan.fa 

cat chr19_BXD005+BXD006len500kb.pan.fa  /home/flaviav/data/BXD/UCSC_mm10_chr19_only.fa > chr19_BXD005+BXD006len1M.pan+ref.fa 
```
I repeated the above commands for DBA2J and the 6 similar strains.

####  4) I used pggb (chr19)

*flaviav@penguin2:/home/flaviav/git/pggb/chr19_DBA2Jstrains_1M_p98*

```shell
pggb -t 30 -i  chr19_DBA2J.pan+ref.fa.gz -Y "#" -p 98 -s 50000 -l 200000 -n 10 -k 29 -B 10000000 -w 10000000 -G 5000 -v -o chr19_DBA2Jstrains_1M_p98

```
#### 5) Variant calling with vg on the pangenome of chr19

```shell
vg deconstruct -e -a -p REF#chr19 -A BXD006 -A BXD008 -A BXD038 -A BXD0147 -A BXD198 -A BXD213 -A DBA2J chr19_DBA2J.pan+ref.fa.gz.2807bd8.2
ff309f.891e76b.smooth.gfa > chr19_DBA2Jstrains.pan+ref.vcf
```

#### 5) I Normalized and sort the VCF file

```shell

bcftools norm --multiallelics -both -f /home/flaviav/data/BXD/UCSC_mm10_chr19_only.fa -o chr19_DBA2Jstrains.pan+ref.norm.decom.vcf chr19_DBA2Jstrains.pan+ref.v
cf
bcftools sort chr19_DBA2Jstrains.pan+ref.norm.decom.vcf > chr19_DBA2Jstrains.pan+ref.norm.decom.sort.vcf
vim chr19_DBA2Jstrains.pan+ref.norm.decom.sort.vcf  #%s/REF#chr19/chr19/g---> for change ID of Reference, for the statistics between genomic and pangenomic VCF, the IDs of reference should be the same.
```
#### 5) I extracted only one sample from the pangenomic VCF, to compare better with the genomic VCF that contains only this sample. I removed sites with 0/0 and positions that have Ns as ALT and REF alleles.
I checked the numbers of positions that have Ns with this: [numbersofNsfromVCFfiles.R](https://github.com/Flavia95/pangmouse/blob/main/script/checkvcfNs.R)

```
bcftools view -s DBA2J  chr19_DBA2Jstrains.pan+ref.norm.decom.sort.vcf > chr19_DBA2Jstrains.pan+ref.norm.decom.sort.onlyDBA2J.vcf
bcftools view -c1 chr19_DBA2Jstrains.pan+ref.norm.decom.sort.onlyDBA2J.vcf > chr19_DBA2Jstrains.pan+ref.norm.decom.sort.onlyDBA2J.filter.vcf
awk '$5  $4 !~ /N/ {print $ALL}' chr19_DBA2Jstrains.pan+ref.norm.decom.sort.onlyDBA2J.filter.vcf > chr19_DBA2Jstrains.pan+ref.norm.decom.sort.onlyDBA2J.filter.withoutN.vcf

```
#### 5) Stats on VCF files

[statspnagvsgeno.sh](https://github.com/Flavia95/pangmouse/blob/main/script/statspangvsgeno.sh)

[allstatistics.R](https://github.com/Flavia95/pangmouse/blob/main/script/allstatistics.R)

---------------------------------------------------------------------------------------------------------------------------------------------------
Other stats:
```shell
bcftools stats BXD005_BXD006.norm.uniq.decomp.vcf >BXD005_BXD006.stats
plot-vcfstats -p outdir BXD005_BXD006.stats  #result-->outdir
vcflib vcflength BXD005_BXD006.norm.uniq.decomp.vcf | vcfbreakmulti | vcf2tsv | cut -f 16 | cut -f 1- >distribution_indelsvcflib.tsv
bcftools view -v other /home/flaviav/data/BXD/BXD005+BXD006_1M/variantcalling/BXD005_BXD006.norm.uniq.decomp.vcf > only_other.vcf
bcftools view -v snps /home/flaviav/data/BXD/BXD005+BXD006_1M/variantcalling/BXD005_BXD006.norm.uniq.decomp.vcf > only_snps.vcf
bcftools view -v indels /home/flaviav/data/BXD/BXD005+BXD006_1M/variantcalling/BXD005_BXD006.norm.uniq.decomp.vcf > only_indels.vcf
bcftools view -v mnps /home/flaviav/data/BXD/BXD005+BXD006_1M/variantcalling/BXD005_BXD006.norm.uniq.decomp.vcf > only_mnp.vcf #variants_classification.tsv (I joined all files and I added column with type of variant)
 ```
**Summary of the statistics**:

1. From variants.txt [script.R](https://github.com/Flavia95/Rplots/blob/main/script/piechart.R):
 [Numbers of variants.png](https://github.com/Flavia95/pangmouse/blob/main/img/distributiononpangenome.png)

2. From distribution_indelsvcflib.tsv [script.R](https://github.com/Flavia95/Rplots/blob/main/script/distributionindels.R):
[Distribution of indels.png](https://github.com/Flavia95/pangmouse/blob/main/img/Distributionofindels.png)

3. From variants_classification.tsv [script.R](https://github.com/Flavia95/Rplots/blob/main/script/variantsalongchromosome.R):
[Distribution variants on chr19.png](https://github.com/Flavia95/pangmouse/blob/main/img/Distributionofvariants.png) 

4. From outdir directory [Substitution types.png](https://github.com/Flavia95/pangmouse/blob/main/img/ts_tv.png)

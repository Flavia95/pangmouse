### Pangenome with only contigs > 1M
I started with two strains.

*flaviav@penguin2:/home/flaviav/data/BXD/BXD005+BXD006_1M*

A. I extract only the contigs greater than 1 million
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
####  4) I used pggb (chr19)

*flaviav@penguin2:/home/flaviav/data/BXD/BXD005+BXD006_1M/outpggb*

```shell
pggb -i /home/flaviav/data/BXD/BXD005+BXD006_1M/whole/chr19_BXD005+BXD006len1M.pan+ref.fa -W -s 50000 -l 150000 -p 90 -w 500000 -j 1200 -e 12000 -n 10 -t 30 -v
 -Y "#" -k 191 -B 2100000 -I 0.95 -R 0.3 --poa-params 1,9,16,2,41,1 -o pggb_BXD005+BXD006_1M

```
#### 5) Variant calling with vg on the pangenome of chr19

*flaviav@penguin2:/home/flaviav/data/BXD/BXD005+BXD006_1M/variantcalling*
```shell
vg deconstruct -e -a -p REF#chr19 -A BXD006 -A BXD005 /home/flaviav/data/BXD/BXD005+BXD006_1M/outpggb/pggb_BXD005+BXD006_1M/chr19_BXD005+BXD006len1Mkb.pan+ref.f
a.pggb-W-s50000-l150000-p90-n10-a0-K16.seqwish-k191-B2100000.smooth-w500000-j1200-e12000-I0.95-p1_9_16_2_41_1.gfa > /home/flaviav/data/BXD/BXD005+BXD006_1M/variantcalling/BXD005_BXD006.vcf
```

#### 5) I Normalized and decomposed the VCF file

*flaviav@penguin2:/home/flaviav/data/BXD/BXD005+BXD006_1M/variantcalling*

```shell

vt index BXD005_BXD006.vcf
vt normalize BXD005_BXD006.vcf.gz -r /home/flaviav/data/BXD/UCSC_mm10_chr19_only.fa | vt uniq - -o BXD005_BXD006.norm.uniq.vcf
vt decompose BXD005_BXD006.norm.uniq.vcf -o BXD005_BXD006.norm.uniq.decomp.vcf
vt peek BXD005_BXD006.norm.uniq.decomp.vcf
```
variants.txt:

Variants          | Number       
--------------| -------------  
SNP        | 104,550  
MNP    | 	6,022
INDEL  | 35,154

#### 5) Stats on VCF file (pangenome chr19)

*flaviav@penguin2:/home/flaviav/data/BXD/BXD005+BXD006_1M/variantcalling/stats*
```shell
bcftools stats BXD005_BXD006.norm.uniq.decomp.vcf >BXD005_BXD006.stats
plot-vcfstats -p outdir BXD005_BXD006.stats  #result-->outdir
vcflib vcflength BXD005_BXD006.norm.uniq.decomp.vcf | vcfbreakmulti | vcf2tsv | cut -f 16 | cut -f 1- >distribution_indelsvcflib.tsv
bcftools view -v other /home/flaviav/data/BXD/BXD005+BXD006_1M/variantcalling/BXD005_BXD006.norm.uniq.decomp.vcf > only_other.vcf
bcftools view -v snps /home/flaviav/data/BXD/BXD005+BXD006_1M/variantcalling/BXD005_BXD006.norm.uniq.decomp.vcf > only_snps.vcf
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

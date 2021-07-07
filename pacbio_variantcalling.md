## Variant calling (Structural Variants) on pacbio data

*flaviav@penguin2:/home/flaviav/data/pacbio*

3_C01--> BXD32; 
5_E01--> B6; 
6_F01--> D2;
7_G01--> BXD1

#### 1. Mapping unaligned BAM files vs reference genome (UCSC_mm10.fa) using [pbmm2](https://github.com/PacificBiosciences/pbmm2)
```shell
pbmm2 align /home/flaviav/data/pacbio/UCSC_mm10.mmi 3_C01/m64247_210429_111231.ccs.bam --sort -j 4 -J 2 3_C01/m64247_210429_111231.ccsVSref.bam
pbmm2 align /home/flaviav/data/pacbio/UCSC_mm10.mmi 5_E01/m64247_210421_013840.ccs.bam --sort -j 4 -J 2 5_E01/m64247_210421_013840.ccsVSref.bam
pbmm2 align /home/flaviav/data/pacbio/UCSC_mm10.mmi 6_F01/m64247_210422_080720.ccs.bam --sort -j 4 -J 2 6_F01/m64247_210422_080720.ccsVSref.bam
pbmm2 align /home/flaviav/data/pacbio/UCSC_mm10.mmi 7_G01/m64247_210423_141129.ccs.bam --sort -j 4 -J 2 7_G01/m64247_210423_141129.ccsVSref.bam
```
#### 2. Adding MD tag, the variant caller "Sniffles" requires alignments to contain the ‘MD tag’ in our BAM file. This is a condensed representation of the alignment of a read to the reference, and is similar to a CIGAR string. Using [samtools calmd](http://www.htslib.org/doc/samtools-calmd.html)
```shell
samtools calmd 3_C01/m64247_210429_111231.ccsVSref.bam -b > 3_C01/m64247_210429_111231.ccsVSref.calmd.bam
samtools calmd 5_E01/m64247_210421_013840.ccsVSref.calmd.bam -b > 5_E01/m64247_210421_013840.ccsVSref.calmd.bam
samtools calmd 6_F01/m64247_210422_080720.ccsVSref.calmd.bam -b > 6_F01/m64247_210422_080720.ccsVSref.calmd.bam
samtools calmd 7_G01/m64247_210423_141129.ccsVSref.calmd.bam -b > 7_G01/m64247_210423_141129.ccsVSref.calmd.bam
```
#### 3. Indexing BAM files, using [samtools index](http://www.htslib.org/doc/samtools-index.html)
```shell
samtools index 3_C01/m64247_210429_111231.ccsVSref.calmd.bam
samtools index 5_E01/m64247_210421_013840.ccsVSref.calmd.bam
samtools index 6_F01/m64247_210422_080720.ccsVSref.calmd.bam
samtools index 7_G01/m64247_210423_141129.ccsVSref.calmd.bam
```
#### 4. [Sniffles](https://github.com/fritzsedlazeck/Sniffles) variant calling
```shell
/home/flaviav/tools/Sniffles-master/bin/sniffles-core-1.0.12/sniffles --genotype -m 3_C01/m64247_210429_111231.ccsVSref.calmd.bam -v 3_C01/m64247_210429_111231.SV.vcf
/home/flaviav/tools/Sniffles-master/bin/sniffles-core-1.0.12/sniffles --genotype -m 5_E01/m64247_210421_013840.ccsVSref.calmd.bam -v /home/flaviav/data/pacbio/5_E01/m64247_210421_013840.SV.vcf
/home/flaviav/tools/Sniffles-master/bin/sniffles-core-1.0.12/sniffles --genotype -m 6_F01/m64247_210422_080720.ccsVSref.calmd.bam -v /home/flaviav/data/pacbio/6_F01/m64247_210422_080720.SV.vcf
/home/flaviav/tools/Sniffles-master/bin/sniffles-core-1.0.12/sniffles --genotype -m 7_G01/m64247_210423_141129.ccsVSref.bam -v /home/flaviav/data/pacbio/7_G01/m64247_210423_141129.SV.vcf
```
Sniffles has a default setting called ‘read support’ which requires 10 reads to support a possible SV for it to be accepted as genuine. 
Reducing this number allows more SVs to be discovered, but may also cause some false positives. **We could try this, for example minimum Support: 5**.

#### 5. Sorting VCF files (crucial step), using [vcf-sort](http://vcftools.sourceforge.net/man_latest.html)
```shell
vcf-sort 3_C01/m64247_210429_111231.SV.vcf > 3_C01/m64247_210429_111231.SV.sort.vcf
vcf-sort 5_E01/m64247_210421_013840.SV.vcf > 5_E01/m64247_210421_013840.SV.sort.vcf
vcf-sort 6_F01/m64247_210422_080720.SV.vcf > 6_F01/m64247_210422_080720.SV.sort.vcf
vcf-sort 7_G01/m64247_210423_141129.SV.vcf > 7_G01/m64247_210423_141129.SV.sort.vcf
```
#### 6. Extracting only chr19 from VCF files keeping the header, using [vcftools](http://vcftools.sourceforge.net/man_latest.html)

vcftools --vcf m64247_210429_111231.SV.sort.vcf --chr chr19 --out c --recode
```shell
vcftools --vcf 3_C01/m64247_210429_111231.SV.sort.vcf --chr chr19 --out 3_C01/m64247_210429_111231.BXD032.SV.chr19 --recode
vcftools --vcf 5_E01/m64247_210421_013840.SV.sort.vcf --chr chr19 --out 5_E01/m64247_210421_013840.BXD6.SV.chr19 --recode
vcftools --vcf 6_F01/m64247_210422_080720.SV.sort.vcf --chr chr19 --out 6_F01/m64247_210422_080720.DBA2J.SV.chr19 --recode
vcftools --vcf 7_G01/m64247_210423_141129.SV.sort.vcf --chr chr19 --out 7_G01/m64247_210423_141129.BXD001.SV.chr19 --recode
```
- Adjusting names, bgzip and tabix VCF files.

Final VCF paths:

*flaviav@penguin2:/home/flaviav/data/pacbio/3_C01/m64247_210429_111231.BXD032.SV.chr19.vcf.gz*

*flaviav@penguin2:/home/flaviav/data/pacbio/5_E01/m64247_210421_013840.BXD6.SV.chr19.vcf.gz*

*flaviav@penguin2:/home/flaviav/data/pacbio/6_F01/m64247_210422_080720.DBA2J.SV.chr19.vcf.g*

*flaviav@penguin2:/home/flaviav/data/pacbio/7_G01/m64247_210423_141129.BXD001.SV.chr19.vcf.gz*

#### 7. Merging all samples in one multi-sample VCF file, using [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR/wiki)

I tried to merge all the VCFs with the standard *bcftools merge*; however, since the breakpoints of the structural variants are only estimates, the position of the same structural variant in different samples would be slightly different. For this thing, I used the SURVIVOR merge, which has a distance setting specifying the maximal difference between pairs of breakpoints (begin1 vs begin2, end1 vs end2). The tool doesn't use a percentage of overlap, but it is also not requiring exact breakpoints. Survivor not uses the percentage of overlap, for the author it is not correct that large SVs reach that minimum percentage easily; for this it uses only a distance.
 
The program uses a distance of 1000 bp. **We could try to use 500bp**.
```shell
ls *vcf > mylist
./SURVIVOR merge samplelist.txt 1000 1 1 1 0 30 merged.vcf
merging entries: 32
merging entries: 717
merging entries: 610
merging entries: 487
merging entries: 355
merging entries: 512
```
*flaviav@penguin2:/home/flaviav/data/pacbio/merged_allstrains.sort.vcf.gz*

#### 8. VCF file (merged_allstrains.sort.vcf.gz) annoteted, with [AnnotSV](https://github.com/lgmgeo/AnnotSV/blob/master/commandLineOptions.txt)
 
./AnnotSV -SVinputFile merged_allstrains.sort.vcf.gz -outputDir /home/flaviav/data/pacbio -genomeBuild mm10 -SVminSize 30

*flaviav@penguin2:/home/flaviav/data/pacbio/merged_allstrains.sort.annotated.tsv*

#### 9. Comparison with VCF files obtained by LongRanger (10X technology)--> in progress
This is an important step, because is not necessary look only the position for both VCF files, but if in the start1 and end1 for the first VCF file there is the SV in the other VCF.
- I'm merging 4 VCF files obtained by LongRanger with Survivor.
- I will use this https://github.com/zeeev/mergeSVcallers or SURVIVOR genComp to the comparison between pacbio VCF and 10X VCF files.

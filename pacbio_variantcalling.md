Variant calling (Structural Variants) on pacbio data

*flaviav@penguin2:/home/flaviav/data/pacbio*
3_C01 - BXD32; 5_E01 - B6; 6_F01 - D2; 7_G01 - BXD1

### 1. Mapped unaligned BAM files vs reference genome 

```shell
pbmm2 align /home/flaviav/data/pacbio/UCSC_mm10.mmi m64247_210429_111231.ccs.bam --sort -j 4 -J 2 3_C01/m64247_210429_111231.ccsVSref.bam
pbmm2 align /home/flaviav/data/pacbio/UCSC_mm10.mmi m64247_210421_013840.ccs.bam --sort -j 4 -J 2 5_E01/m64247_210421_013840.ccsVSref.bam
pbmm2 align /home/flaviav/data/pacbio/UCSC_mm10.mmi m64247_210422_080720.ccs.bam --sort -j 4 -J 2 6_F01/m64247_210422_080720.ccsVSref.bam
pbmm2 align /home/flaviav/data/pacbio/UCSC_mm10.mmi m64247_210423_141129.ccs.bam --sort -j 4 -J 2 7_G01/m64247_210423_141129.ccsVSref.bam
```
### 2. Added MD tag, sniffles requires alignments to contain the ‘MD tag’ in our BAM file. This is a condensed representation of the alignment of a read to the reference, and is similar to a CIGAR string

```shell
samtools calmd 3_C01/m64247_210429_111231.ccsVSref.bam -b > 3_C01/m64247_210429_111231.ccsVSref.calmd.bam
samtools calmd 5_E01/m64247_210421_013840.ccsVSref.calmd.bam -b > 5_E01/m64247_210421_013840.ccsVSref.calmd.bam
samtools calmd 6_F01/m64247_210422_080720.ccsVSref.calmd.bam -b > 6_F01/m64247_210422_080720.ccsVSref.calmd.bam
samtools calmd 7_G01/m64247_210423_141129.ccsVSref.calmd.bam -b > 7_G01/m64247_210423_141129.ccsVSref.calmd.bam
```
### 3. Index BAM files
```shell
samtools index 3_C01/m64247_210429_111231.ccsVSref.calmd.bam
samtools index 5_E01/m64247_210421_013840.ccsVSref.calmd.bam
samtools index 6_F01/m64247_210422_080720.ccsVSref.calmd.bam
samtools index 7_G01/m64247_210423_141129.ccsVSref.calmd.bam
```
### 4. Sniffles variant calling
```shell
/home/flaviav/tools/Sniffles-master/bin/sniffles-core-1.0.12/sniffles --genotype -m 3_C01/m64247_210429_111231.ccsVSref.calmd.bam -v 3_C01/m64247_210429_111231.SV.vcf

/home/flaviav/tools/Sniffles-master/bin/sniffles-core-1.0.12/sniffles --genotype -m 5_E01/m64247_210421_013840.ccsVSref.calmd.bam -v /home/flaviav/data/pacbio/5_E01/m64247_210421_013840.SV.vcf

/home/flaviav/tools/Sniffles-master/bin/sniffles-core-1.0.12/sniffles --genotype -m 6_F01/m64247_210422_080720.ccsVSref.calmd.bam -v /home/flaviav/data/pacbio/6_F01/m64247_210422_080720.SV.vcf

/home/flaviav/tools/Sniffles-master/bin/sniffles-core-1.0.12/sniffles --genotype -m 7_G01/m64247_210423_141129.ccsVSref.bam -v /home/flaviav/data/pacbio/7_G01/m64247_210423_141129.SV.vcf
```
Sniffles has a default setting called ‘read support’ which requires 10 reads to support a possible SV for it to be accepted as genuine. 
Reducing this number allows more SVs to be discovered, but may also cause some false positives. We could try to re-set.

### 5.Sort VCF files (crucial step)
```shell

/genotype -m 3_C01/m64247_210429_111231.ccsVSref.calmd.bam -v 3_C01/m64247_210429_111231.SV.vcf

/home/flaviav/tools/Sniffles-master/bin/sniffles-core-1.0.12/sniffles --genotype -m 5_E01/m64247_210421_013840.ccsVSref.calmd.bam -v /home/flaviav/data/pacbio/5_E01/m64247_210421_013840.SV.vcf

/home/flaviav/tools/Sniffles-master/bin/sniffles-core-1.0.12/sniffles --genotype -m 6_F01/m64247_210422_080720.ccsVSref.calmd.bam -v /home/flaviav/data/pacbio/6_F01/m64247_210422_080720.SV.vcf


vcf-sort m64247_210423_141129.SV.vcf > m64247_210423_141129.SV.sort.vcf
```



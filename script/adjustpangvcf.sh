#!/bin/bash
VCF=$1
REF="/home/flaviav/data/BXD/148strains/UCSC_mm10_changeid.fa"

#Normalize and remove N alleles
bcftools norm -c s -m - -f ${REF} -o chr19.pan+ref.norm.vcf ${VCF}
awk '$5  $4 !~ /N/ {print $ALL}' chr19.pan+ref.norm.vcf > chr19.pan+ref.filter.vcf && bcftools sort chr19.pan+ref.filter.vcf >chr19.pan+ref.filter.sort.vcf

#Rename samples id and REF
cat chr19.pan+ref.filter.sort.vcf | sed 's/^REF#//g' >chr19.pan+ref.tmp.vcf && bcftools query -l chr19.pan+ref.tmp.vcf | sed 's/#//g' > renamedsamples.txt && bcftools reheader -s renamedsamples.txt chr19.pan+ref.tmp.vcf -o chr19.pan+ref.renamed.vcf && bgzip chr19.pan+ref.renamed.vcf && tabix chr19.pan+ref.renamed.vcf.gz &&
rm chr19.pan+ref.tmp.vcf

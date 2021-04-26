echo "Genomic and pangenomic files: $1 $2"
    
#1.Stats for variable sites

vcf-compare $1 $2 | grep ^VN | cut -f 2- | column -t > vcfcomparepos.txt
vcf-compare $1 $2 | grep ^SN | cut -f 2- | column -t>vcfcomparealtref.txt
vcftoolz compare $1 $2 > vcftoolzcompare.txt

#2. Stats for variants

bcftools stats $1 $2 > bcftools.txt
cat bcftools.txt | grep "^SN\|^# SN" | grep -v "SN," > overlapvariants.tsv

#Indels
vcftools --gzvcf $1 --out geno --hist-indel-len
vcftools --gzvcf $2 --out pang --hist-indel-len

#Positive and negative
#gatk Concordance -R /home/flaviav/data/BXD/variantcalling_genomic/UCSC_mm10_chr19_only.fa  --truth $1 --evaluation $2 --summary summarygatkconcordance.tsv

####bonus_stats

vcftools --gzvcf $1 --freq2 --out freq_geno --max-alleles 2

vcftools --gzvcf $1 --het --out het_gen

vcftools --gzvcf $1 --site-quality --out qual_sites_gen


vcftools --gzvcf $2 --site-quality --out qual_sites_pan

vcftools --gzvcf $2 --freq2 --out freq_pan --max-alleles 2

vcftools --gzvcf $2 --het --out het_pan


#Rscript allstatistics.R

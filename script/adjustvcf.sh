echo "Genomic VCF: $1"

#Normalized Genomic VCF
bcftools norm --multiallelics -both -f /home/davida/UCSC_mm10.fa genomic.vcf -o Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.recode.norm.de.vcf $1

#Reheader Genomic VCF
bcftools reheader -s listrename.txt Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.recode.norm.de.vcf -o Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.reheader.vcf.gz

#Remove spanning deletions
zgrep -v "*" Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.reheader.vcf.gz | /home/flaviav/tools/bgzip > Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.filter.vcf.gz && /home/flaviav/tools/tabix Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.filter.vcf.gz

#Extract only chr19
bcftools view Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.filter.vcf.gz --regions chr19 > Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.filter.chr19.vcf

#Remove samples which are not assembled
vcftools --remove removeind.txt --vcf Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.filter.chr19.vcf --recode --out chr19_removsamples/Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.filter.chr19.removeind.vcf

#Remove sites ref/ref
bcftools view -c1 Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.filter.chr19.removeind.vcf > chr19_removsamples/Merged_gvcf_files_autosomes_Feb_2021_recalibrated_all_variants_PASSED.filter.chr19.removeind_sites.vcf

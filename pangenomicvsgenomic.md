### Comparison between the VCF obtained from genomic model and VCF obtained from pangenomic models.
The VCF files contain three mouse samples: DBA_2J_0334_three_lanes, DBA2_J and DBA_2J.


**A) I use little region---> 3,117,753 - 3,336,924**

#### 1) I apply [vcflift.py](https://github.com/Flavia95/VGpop/blob/master/doc/vcflift.md)
 
It takes as input a chunk of a vcf file, a chunk of the reference sequence in fasta that covers the region in the vcf and adjust the chromosome positions to match the positions in the fasta. The vcf should have been generated from the same reference sequence from which the fasta is extracted.

*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/little_region*


```shell
bedtools getfasta -fi UCSC_mm10.fa -bed region.zerobased.bed > UCSC_mm10.zerobased.fa
```
*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/little_region/snps*
```shell
tabix -h /home/flaviav/testgenomic_pangenomic/DBA_2J_consensus_site_copy_hom_only_chr19_only.vcf.gz chr19:3117753-3336924 > DBA2J_extractregion.vcf

vcflift.py -fasta  /home/flaviav/testgenomic_pangenomic/little_region/UCSC_mm10.zerobased.fa -start 3117752 -vcf DBA2J_extractregion.vcf > DBA2J_consensus_sites_chr19_exractregion_vcflift.vcf 
```

#### 2) I apply vg tool for extract little region from original VCF and call SNPs

*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/little_region/snps/vg_gfautilcallsnps*

```shell

-vg construct -r /home/flaviav/testgenomic_pangenomic/UCSC_mm10_onlychr19.fa -v /home/flaviav/testgenomic_pangenomic/DBA_2J_consensus_site_copy_hom_only_chr19_only.vcf.g -R chr19::3117753-3336924 > DBA2J_chr19_3117753-3336924.vg

-vg view chr19_3117753-3336924.indels.vg> DBA2J_chr19_3117753-3336924.gfa

-vg index DBA2J_chr19_3117753-3336924.indels.vg -x DBA2J_chr19_3117753-3336924.xg

-vg deconstruct -p "chr19" DBA2J_chr19_3117753-3336924.xg > /home/flaviav/testgenomic_pangenomic/little_region/snps/DBA2J_chr19_3117753-3336924.vgcall.vcf

```

####  3) I apply gfautil for call variants

gfautil doesn't work on GFA with one path, for this I add random paths with odgi.

*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/little_region/snps/vg_gfautilcallsnps*

```shell

- odgi build -g DBA2J_chr19_3117753-3336924.gfa -o DBA2J_chr19_3117753-3336924.og

- odgi cover -i DBA2J_chr19_3117753-3336924.og -o DBA2J_chr19_3117753-3336924.allpaths.og  -c 1

- odgi view -i  DBA2J_chr19_3117753-3336924.allpaths.og  -g >  DBA2J_chr19_3117753-3336924.allpaths.gfa

gfautil --debug -t 20 -i  DBA2J_chr19_3117753-3336924.allpaths.gfa  gfa2vcf --refs "chr19"> /home/flaviav/testgenomic_pangenomic/little_region/snps/DBA2J_chr19_3117753-3336924.gfautil.vcf
```

####  4) I compare the output obtained by vcflift.py with adjust the chromosome positions in the VCF file to match the positions in the fasta (genomic model) and output obtained by VG and GFAUTIL (pangenomic models).

*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/little_region/snps/compareDBA2J_DBA2vgcall*

```shell
bcftools isec DBA2J_chr19_3117753-3336924.vgcall.sort.vcf.gz DBA2J_consensus_sites_chr19_exractregion_vcflift.sort.vcf.gz -p compareDBA2J_DBA2vgcall
```

Models           | VCFfiles |SNPs        
--------------| ----|---------  
Genomic model (freebayes)      | DBA2J_consensus_sites_chr19_exractregion_vcflift.sort.vcf.gz|  9
Pangenomic models (vg) | DBA2J_chr19_3117753-3336924.vgcall.sort.vcf.gz| 7
Pangenomic models (gfautil) | DBA2J_chr19_3117753-3336924.gfautil.vcf | 7

The number and positions of the variants is the same for vg, gfautil and the genomic method. 

The genomic method calls two more variants


#### 5) I repeat the same steps for the indels. 

Gfautil fails to call variants even if I add paths with odgi cover (flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/little_region/indels/vg_gfautilcall).

~~Error: (Parsing GFA from chr19_3117753-7117753.indels.allpaths.gfa thread 'main' panicked at 'GFA must contain at least two paths'.)~~ 


*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/little_region/indels*

Models           | VCFfiles |INDELs        
--------------| ----|---------  
Genomic model  (freebayes)  | DBA2J_consensus_indels_chr19_exractregion_vcflift.sort.vcf.gz|  6
Pangenomic models (vg) | DBA2J_chr19_3117753-3336924.indels.vgcall.sort.vcf.gz| 6

The number and positions of the variants is the same for vg and the genomic method. 

-------------------------------------------------------------------------------------------------------------------------------------------------------------
**B) I use a region of 4 MB --> 3,117,752 - 7,117,753**

*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/big_region*

#### 1) I apply vcflift.py 
#### 2) I use vg and gfautil for call variants
#### 3) I compare the output obtained by vcflift.py (genomic model) and output obtained by vg and gfautil (pangenomic model)

- For SNPs

*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/big_region/snps*

Models           | VCFfiles |SNPs        
--------------| ----|---------  
Genomic model (freebayes) | DBA2J_consensus_sites_chr19_exractregion_vcflift.sort.vcf.gz |  3816
Pangenomic models (vg) | DBA2J_chr19_3117753-7117753.vgcall.sort.vcf.gz| 3815
Pangenomic models (gfautil) | DBA2J_chr19_3117753-7117753.gfautil.vcf.gz | 3815

The number and positions of the variants is the same for vg, gfautil and the genomic method. 

The genomic method calls one more variant.

To see which variant is not shared between the genomic method and the pangenomic methods:

*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/big_region/snps/compareDBA2Jvcfliftvsvgcall/uniqvariants*

```shell
bcftools isec -p uniqvariants -n-1 -c all DBA2J_consensus_sites_chr19_exractregion_vcflift.sort.vcf.gz  DBA2J_chr19_3117753-7117753.vgcall.sort.vcf.gz
```

- For INDELs

*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/big_region/indels*

Models           | VCFfiles |INDELs      
--------------| ----|---------  
Genomic model (freebayes)     | DBA2J_consensus_indels_chr19_exractregion_vcflift.sort.vcf.gz |  929
Pangenomic models (vg) | DBA2J_3117753-7117753.indels.vgcall.sort.vcf.gz | 927

flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/big_region/indels/comparetwovcfindels/uniqvariantsindels
To see which variant is not shared between the genomic method and the pangenomic method.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------

**C) ALL CHROMOSOME 19 ---> for this I use only gfautil because vg doesn't work**

*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/test_allchr*

-SNPs

1. I USE BCFTOOLS FOR EXTRACT FASTA SEQUENCES (three samples) FROM VCF file.

```shell

bcftools consensus -i -s DBA_2J_0334_three_lanes -f UCSC_mm10_onlychr19.fa  DBA_2J_consensus_sites_copy_hom_only_chr19_only.vcf.gz  > DBA_2J_0334_three_lanes.fa

bcftools consensus -i -s DBA2_J -f UCSC_mm10_onlychr19.fa  DBA_2J_consensus_sites_copy_hom_only_chr19_only.vcf.gz  > DBA2_J.fa

bcftools consensus -i -s  DBA_2J -f UCSC_mm10_onlychr19.fa  DBA_2J_consensus_sites_copy_hom_only_chr19_only.vcf.gz  >  DBA_2J.fa

```

2. I CHANGE SEQUENCES IDs

```shell
awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"DBA2_J#"$1} else {print $0}}' DBA2_J.fa | grep -o '^\S*' > DBA2_J_changeid.fa

awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"DBA_2J#"$1} else {print $0}}' DBA_2J.fa | grep -o '^\S*' > DBA_2J_changeid.fa

awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"DBA_2J0334#"$1} else {print $0}}' DBA_2J_0334_three_lanes.fa | grep -o '^\S*' > DBA_2J_0334_changeid.fa
```

3. I JOIN 4 SEQUENCES
```shell
cat  DBA2_J_changeid.fa DBA_2J_changeid.fa DBA_2J_0334_changeid.fa UCSC_mm10_onlychr19.changeid.fa > DBA2_J+DBA_2J+DBA_2J_0334+UCSC_chr19.sites.fa  
```
*flaviav@penguin2:/home/flaviav/testgenomic_pangenomic/test_allchr/outDBA2J_treelines*

4. I USE PGGB

```shell

./pggb -i DBA2_J+DBA_2J+DBA_2J_0334+UCSC_chr19.sites.fa -W -a 0 -s 500 -l 1000 -p 95 -w 500000 -j 12000 -e 12000 -n 15 -t 20 -v -Y "#" -k 27 -B 25000000 -I 0.7 -R 0.2 -C 10000 --poa-params 1,9,16,2,41,1 -o outDBA2J_treelines
```

5.  I USE GFAUTIL
 
```shell
gfautil --debug -t 20 -i DBA2_J+DBA_2J+DBA_2J_0334+UCSC_chr19.sites.fa.pggb-W-s500-l1000-p95-n15-a0-K16.seqwish-k27-B25000000.smooth-w500000-j12000-e12000-I0.7-p1_9_16_2_41_1.gfa gfa2vcf --refs "REF#chr19"> DBA2_J+DBA_2J+DBA_2J_0334+refchr19.vcf
```

On this output there are 100 SNPs.


5. I remove paths consensus and redo gfautil.

```shell
grep Consensus DBA2_J+DBA_2J+DBA_2J_0334+UCSC_chr19.sites.fa.pggb-W-s500-l1000-p95-n15-a0-K16.seqwish-k27-B25000000.smooth-w500000-j12000-e12000-I0.7-p1_9_16_2_41_1.gfa -v > DBA2_J+DBA_2J+DBA_2J_0334+UCSC_chr19.sites_without.consensus.gfa


gfautil --debug -t 20 -i DBA2_J+DBA_2J+DBA_2J_0334+UCSC_chr19.sites_without.consensus.gfa gfa2vcf --refs "REF#chr19"> DBA2_J+DBA_2J+DBA_2J_0334+UCSC_chr19.sites_without.consensus.vcf

```
 On this output there are 15 variants.



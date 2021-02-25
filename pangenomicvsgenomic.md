
1. BCFTOOLS FOR EXTRACT FASTA FILES FROM VCF

```shell

bcftools consensus -i -s DBA_2J_0334_three_lanes -f UCSC_mm10_onlychr19.fa  DBA_2J_consensus_sites_copy_hom_only_chr19_only.vcf.gz  > DBA_2J_0334_three_lanes.fa

bcftools consensus -i -s DBA2_J -f UCSC_mm10_onlychr19.fa  DBA_2J_consensus_sites_copy_hom_only_chr19_only.vcf.gz  > DBA2_J.fa

bcftools consensus -i -s  DBA_2J -f UCSC_mm10_onlychr19.fa  DBA_2J_consensus_sites_copy_hom_only_chr19_only.vcf.gz  >  DBA_2J.fa

```


2. CHANGE SEQUENCES ID

```shell
awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"DBA2_J#"$1} else {print $0}}' DBA2_J.fa | grep -o '^\S*' > DBA2_J_changeid.fa

awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"DBA_2J#"$1} else {print $0}}' DBA_2J.fa | grep -o '^\S*' > DBA_2J_changeid.fa

awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"DBA_2J0334#"$1} else {print $0}}' DBA_2J_0334_three_lanes.fa | grep -o '^\S*' > DBA_2J_0334_changeid.fa
```

3. JOIN 4 SEQUENCES
```shell
cat  DBA2_J_changeid.fa DBA_2J_changeid.fa DBA_2J_0334_changeid.fa UCSC_mm10_onlychr19.changeid.fa > DBA2_J+DBA_2J+DBA_2J_0334+UCSC_chr19.sites.fa  
```

4. PGGB

```shell

./pggb -i DBA2_J+DBA_2J+DBA_2J_0334+UCSC_chr19.sites.fa -W -a 0 -s 500 -l 1000 -p 95 -w 500000 -j 12000 -e 12000 -n 15 -t 20 -v -Y "#" -k 27 -B 25000000 -I 0.7 -R 0.2 -C 10000 --poa-params 1,9,16,2,41,1 -o outDBA2J_treelines
```

5. Remove paths consensus and gfautil

```shell
awk '!/Consensus/' DBA2_J+DBA_2J+DBA_2J_0334+UCSC_chr19.sites.fa.pggb-W-s500-l1000-p95-n15-a0-K16.seqwish-k27-B25000000.smooth-w500000-j12000-e12000-I0.7-p1_9_16_2_41_1.gfa > DBA2_J+DBA_2J+DBA_2J_033.withoutconses.gfa

gfautil --debug -t 20 -i DBA2_J+DBA_2J+DBA_2J_033.withoutconses.gfa gfa2vcf --refs "REF#chr19"> prova.vcf

```

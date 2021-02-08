## Build pangenome for more strains of mouse

I started from assembly obtained by Supernova tool

-  For visualized a compact nameID then in pangenome I changed ID for each sequence in each strain:

```shell
awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"DBA2J#"$1} else {print $0}}' /home/flaviav/data/DBA2J_supernova_3lanes_de_novo.fasta | grep -o '^\S*' > DBA2J_supernova_changeid.fasta
```

- Join fasta files with ID changed

```shell
cat BXD001_supernova_changeid.fasta BXD002_supernova_changeid.fasta BXD005_supernova_changeid.fasta BXD006_supernova_changeid.fasta BXD008_supernova_changeid.fasta BXD009_supernova_changeid.fasta BXD011_supernova_changeid.fasta BXD012_supernova_changeid.fasta C57BL6J_supernova_changeid.fasta DBA2J_supernova_changeid.fasta > C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_supernova_changeid.fa
```

- Mapping all strains against whole reference 

```shell
time edyeet -m -N -p 95 -Y '#' -t 20 -Q <(ls C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_supernova_changeid.fa.gz | sort -V ) UCSC_mm10_changeid.fa > C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012vsref.paf
```

- Extracted IDs sequences of chr19

```shell
echo 19 | while read i; do awk '$6 == "REF#chr'$i'"' C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012vsref.paf; done | cut -f 1 | sort -V > id_C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012vsref.chr19.txt
```

- Ectracted only rows of chr19 from fasta files (on not gzip file)

```shell 
cat id_C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012vsref.chr19.txt | while read line ; do samtools faidx C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_supernova_changeid.fa $line; done > C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_chr19.fa
```
- pggb pipeline

```shell
./pggb -i /home/flaviav/C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_chr19.fa.gz -s 5000 -p 99 -n 10 -Y '#' -w 30000 -s 50000 -I 0.5 -t 16 -o out
```

## In the left of the pangenome there is an unconstructed portion, that seems correspond to the centromeric region, I checked this

1. I extracted all sequences IDS in the fasta files of ten strains (onlychr19):
```shell
zgrep "^>" C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_chr19.fa.gz | cut -d'|' -f3 | cut -d ' ' -f1 > file.txt
```
2. For each sequence ID I added the position that correspondes to the centromer (0-3MB) and I obtained a bed file

3. I extracted for each sequence IDs the sequences that correspondes to the centromer:
```shell
samtools faidx C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_chr19.fa

bedtools getfasta -fi C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_chr19.fa -bed centromericregion.bed > 10strainsonlypos
centro.fa
```
4. I built with pggb the pangenome --> it is a zoom of the left of the previous pangenome and it correpondes to the centromeric region..

### VARIANT CALLING ON THE PANGENOME OBTAINED BY PGGB

## The best command for build pangenome from 39 strains:

```shell
./pggb -i /home/flaviav/data/BXD/39strains_chr19+ref_chr19.fa -a 0 -s 5000 -l 15000 -p 99 -w 30000 -j 5000  -e 5000 -n 6 -t 20 -Y "#" -k 29 -K 16 -I 0.5 -R 0.2 -o pang5000

./pggb -i /home/flaviav/data/BXD/39strains_chr19+ref_chr19.fa -a 0 -s 8000 -p 98 -w 300000 -j 8000 -e 8000 -n 6 -t 20 -Y "#" -S -k 29 -I 0.6 -R 0.2 --resume -o pang8000
```

I'm based on the output obtained by the last command:

- Stastics on GFA obtained by seqwish: 
#length   nodes   edges   paths
257572767 1406518 2521219 19773

- Statistics on GFA obtained by smoothxg:
#length    nodes   edges     paths
223294493  9891933 14023885  99633


- Variant calling

```shell
env TMPDIR=~/data/tmp gfautil --debug -t 20 -i /home/flaviav/data/pang_39strains/pang8000/39strains_chr19+ref_chr19.fa.pggb-E-s8000-l24000-p98-n6-a0-K16-k29-w300000-j8000-e8000-I0.6-R0.2.smooth.gfa gfa2vcf --refs "REF#chr19" > /home/flaviav/data/pang_39strains/call_variants/39strains_s8000_smooth_pangenome.vcf
```

- Normalized the output of gfautil
```shell
env TMPDIR=~/data/tmp vt index 39strains_s8000_smooth_pangenome.vcf.gz

env TMPDIR=~/data/tmp vt normalize -n 39strains_s8000_smooth_pangenome.vcf.gz -r /home/flaviav/data/UCSC_mm10_chr19_only.fa | vt uniq - -o 39strains_s8000_smooth_pangenome.norm.uniq.vcf

env TMPDIR=~/data/tmp vt decompose 39strains_s8000_smooth_pangenome.norm.uniq.vcf -o 39strains_s8000_smooth_pangenome.norm.uniq.decomp.vcf
```

Plot

```R
myd = read.table("/home/flavia/Desktop/table.txt", header=T)
p = myd %>% gather(count, value, values:N) %>% ggplot(aes(value, samples, fill = count)) + geom_bar(stat = "identity", width=0.50) + facet_wrap(.~ species)
p + labs(x="Values") + theme(legend.title = element_blank()) + theme_bw()+ + ylab("Chromosomes") + scale_fill_discrete(name = “Count”, labels = c(“N”, "nucleotides"))

ggsave("/home/flavia/Desktop/Distribution", plot= p, device="png", width = 20, height = 15, units = "cm", dpi = 300)
```

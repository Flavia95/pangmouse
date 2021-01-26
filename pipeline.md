## Build pangenome for more strains of mouse

I started from assembly obtained by Supernova tool

1. For visualized a compact nameID then in pangenome I changed ID for each sequence in each strain:

```shell
awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"DBA2J#"$1} else {print $0}}' /home/flaviav/data/DBA2J_supernova_3lanes_de_novo.fasta | grep -o '^\S*' > DBA2J_supernova_changeid.fasta
```

2. Join fasta files with ID changed

```shell
cat BXD001_supernova_changeid.fasta BXD002_supernova_changeid.fasta BXD005_supernova_changeid.fasta BXD006_supernova_changeid.fasta BXD008_supernova_changeid.fasta BXD009_supernova_changeid.fasta BXD011_supernova_changeid.fasta BXD012_supernova_changeid.fasta C57BL6J_supernova_changeid.fasta DBA2J_supernova_changeid.fasta > C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_supernova_changeid.fa
```

3. Mapping all strains against whole reference 

```shell
time edyeet -m -N -p 95 -Y '#' -t 20 -Q <(ls C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_supernova_changeid.fa.gz | sort -V ) UCSC_mm10_changeid.fa > C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012vsref.paf
```

4. Extracted IDs sequences of chr19

```shell
echo 19 | while read i; do awk '$6 == "REF#chr'$i'"' C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012vsref.paf; done | cut -f 1 | sort -V > id_C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012vsref.chr19.txt
```

5. Ectracted only rows of chr19 from fasta files (on not gzip file)

```shell 
cat id_C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012vsref.chr19.txt | while read line ; do samtools faidx C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_supernova_changeid.fa $line; done > C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012_chr19.fa
```
6. pggb pipeline

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


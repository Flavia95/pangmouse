### Cover centromeric region



#### 1) Remove centromeric regions from the reference of mouse

```shell
samtools faidx UCSC_mm10_changeid.fa
awk 'BEGIN { FS="\t" } { print $1, $2 }' UCSC_mm10_changeid.fa.fai > extractposwithoutcentr.bed
bedtools getfasta -fi UCSC_mm10_changeid.fa -bed extractposwithoutcentr.bed > UCSC_mm10_changeid.withoutcentro.fa
```
I changed sequence IDs:

```shell

sed 's/^>.*read/>read/' DBA2J.fa > DBA2J.changeid.fa
sed 's, ,_,g' -i DBA2J.changeid.fa
```

#### 2) MINIMAP THAT WORKS WELL WITH READS NANOPORE

For each strain:


```shell

env TMPDIR=~/data/tmp ./minimap2 -ax map-ont /home/flaviav/data/nanopore_parental/minimap/UCSC_mm10_changeid.withoutcentro.fa /home/flaviav/data/nanopore_parental/DBA2J.changeid.fa > /home/flaviav/data/nanopore_parental/minimap/DBA2J.all.sam

```

####  3) Split mapped and unmapped


```shell
samtools view -F4 DBA2J.all.sam > DBA2J.all.mapped.sam
samtools view -f4 DBA2J.all.sam > DBA2J.all.unmapped.sam

```

####  4) From mapped reads I extracted only reads of chr19

```shell
awk '$3 ~ /REF#chr19/ { print }' DBA2J.all.mapped.sam > DBA2J.all.mapped.chr19.sam
```

#### 5) I extracted only reads that mapped at left

python script

```shell
python3 parserreads.py -input chr19nanopore/DBA2J.all.mapped.chr19.sam -output chr19nanopore/id_DBA2J.all.mapped.chr19.txt
```

#### 6) I extracted only reads that mapped at left

python script

```shell
python3 parserreads.py -input chr19nanopore/DBA2J.all.mapped.chr19.sam -output chr19nanopore/id_DBA2J.all.mapped.chr19.txt
```

```shell
cat -n chr19nanopore/id_DBA2J.all.mapped.chr19.txt | sort -uk2 | sort -nk1 | cut -f2- > chr19nanopore/id_DBA2J.all.mapped.chr19_unique.txt
```

#### 7) I extracted only sequence IDs from fasta file

```shell
samtools faidx /home/flaviav/data/nanopore_parental/C57BL6J.changeid.fa  $(cat chr19nanopore/id_DBA2J.all.mapped.chr19_unique.txt ) > chr19nanopore/DBA2J.all.mapped.chr19.onlyreadscentromer.fa
```

#### 8) Join all

```shell
cat DBA2J.all.mapped.chr19.changeid.pggb.fa C57BL6J.all.mapped.chr19.changeid.pggb.fa /home/flaviav/data/BXD/assemblies/UCSC_mm10_chr19_only.fa > DBA2J.C57BL6J.ref_onlychr19.fa
```
#### 8) PGGB
There are reads that cover centromer but fragmented..

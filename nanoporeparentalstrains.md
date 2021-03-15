### Nanopore reads of two parental strains

Steps before starting:

```shell
paste - - - - < DBA_2J.fq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > DBA_2J.fa
paste - - - - < C57BL_6J.fq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > C57BL6J.fa
```

I fix sequence IDs:

```shell

sed 's/^>.*read/>read/' DBA2J.fa > DBA2J.changeid.fa
sed 's, ,_,g' -i DBA2J.changeid.fa
```

#### 1) I remove centromeric regions from the reference sequence of mouse

```shell
samtools faidx UCSC_mm10_changeid.fa
awk 'BEGIN { FS="\t" } { print $1, $2 }' UCSC_mm10_changeid.fa.fai > extractposwithoutcentr.bed  #I extract len to each chr
bedtools getfasta -fi UCSC_mm10_changeid.fa -bed extractposwithoutcentr.bed > UCSC_mm10_changeid.withoutcentro.fa
```

#### 2) I map each strain against reference genome without centromeric region

Same for C57BL6J
```shell

minimap2 -ax map-ont UCSC_mm10_changeid.withoutcentro.fa DBA2J.changeid.fa > DBA2J.all.sam

```

####  3) I extract from the sam mapped and unmapped reads

```shell
samtools view -F4 DBA2J.all.sam > DBA2J.all.mapped.sam
samtools view -f4 DBA2J.all.sam > DBA2J.all.unmapped.sam

```

####  4) From the file containing the mapped reads I extract only the chr19 reads

```shell
awk '$3 ~ /REF#chr19/ { print }' DBA2J.all.mapped.sam > DBA2J.all.mapped.chr19.sam
```

#### 5) Take the reads that maps to the left of reference genome 

I consider field 4 of the SAM file--> 1-based leftmost mapping POSition.

If the left position is less than or equal to the length of the reads, I extract the IDs of the reads. 

[parserreads.py](script/parserreads.py)

```shell
python3 parserreads.py -input DBA2J.all.mapped.chr19.sam -output id_DBA2J.all.mapped.chr19.txt
```
```shell
cat -n id_DBA2J.all.mapped.chr19.txt | sort -uk2 | sort -nk1 | cut -f2- > id_DBA2J.all.mapped.chr19_unique.txt
```

#### 6) I extract only sequence IDs from fasta file

```shell
samtools faidx DBA2J.changeid.fa  $(cat id_DBA2J.all.mapped.chr19_unique.txt ) > DBA2J.all.mapped.chr19.onlyreadscentromer.fa
```
I change IDs sequences.

#### 7) Join two parental strains and the reference of chr19

```shell
cat DBA2J.all.mapped.chr19.changeid.pggb.fa C57BL6J.all.mapped.chr19.changeid.pggb.fa UCSC_mm10_chr19_only.fa > DBA2J.C57BL6J.ref_onlychr19.fa
```
#### 8) PGGB
There are reads that cover centromeric regions but I need of an assembly.

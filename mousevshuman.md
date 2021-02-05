# Mapping reference of mousevshuman

Step before starting

- ### I changed IDs of reference sequences:

```shell
awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"Mouse#"$1} else {print $0}}' UCSC_mm10.fa | grep -o '^\S*' > UCSC_mm10.fa_changeid.fa
```
```shell
awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"Homosapiens#"$1} else {print $0}}' GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna | grep -o '^\S*' > GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_changeid.fa
```

- ### I tried lastz

https://www.ensembl.org/Homo_sapiens/Location/Synteny?r=17:63992802-64038237;otherspecies=Mus_musculus

I downloaded these sequences:

PRR29 (ENSG00000224383)	17:63998351-64004305	→	Prr29 (ENSMUSG00000009210)	11:106365472-106377558

```shell
lastz chr17hum.fa chr11mouse.fa --notransition --step=20 --nogapped --format=maf > aln.maf
last-dotplot aln.maf algn.png
```
![algn.png](/img/algn.png)

## Approximate mapping

- wfmash

```shell
wfmash -m -K -k 29 -p 75 -s 10000000 -l 0 -t 16 GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_changeid.fa.gz UCSC_mm10.fa_changeid.fa.gz > human-mouse.-m-k29-p75-K-s10000000-l0.paf
```

```shell
sort -n -k 11 human-mouse.-m-k29-p75-K-s10000000-l0.paf | awk '{ sum += $11 } END { print sum }'
1700012320
```

## Possible use of approximate mapping + lastz

1. Taked a single line from the PAF obtained by approximate mapping (wfmash).

I Used the coordinates and strand orientation to extract the query and target subsequence and then I aligned these with lastz.

Ex: #filename: extractseq_mousefrompaf.bed 

Mouse#chrX	30000000	60000000	-                   

#filename: extractseq_humanfrompaf.bed

Homosapiens#chr1	44910958	73919870	-             

```shell
bedtools getfasta -fi mouse/UCSC_mm10_changeid.fa -bed extractseq_mousefrompaf.bed > UCSC_idfrompaf.fa
bedtools getfasta -fi human/GCA_000001405.27_GRCh38.p12_genomic_changeid.fa -bed extractseq_humanfrompaf.bed > GCA_idfrompaf.fa
```

```shell
lastz GCA_idfrompaf.fa UCSC_idfrompaf.fa --notransition --step=20 --nogapped --format=maf > aln.maf
last-dotplot aln.maf algnChrxandChr1.png
```

![algnChrxandChr1.png](/img/algnChrxandChr1.png)
or with R

```shell
lastz GCA_idfrompaf.fa UCSC_idfrompaf.fa --notransition --step=20 --nogapped --format=rdotplot >  algnChrxandChr1.txt
```
```R
dots = read.table("algnChrxandChr1.txt",header=T,row.names=NULL)
plot(dots,type="l")
```
![sample.png](/img/sample.png)

or in R with my script: 
[dotplot.R](script/dotplot.R)

![alnchrx_chr1.png](/img/alnchrx_chr1.png)

## I tried with more IDS:

From these bed files:

[extractseq_humanfrompaf](test/extractseq_humanfrompaf.bed)

[extractseqmousefrompaf](test/extractseq_mousefrompaf.bed)

```shell
bedtools getfasta -fi mouse/UCSC_mm10_changeid.fa -bed extractseq_mousefrompaf.bed > UCSC_idfrompaf.fa
bedtools getfasta -fi human/GCA_000001405.27_GRCh38.p12_genomic_changeid.fa -bed extractseq_humanfrompaf.bed > GCA_idfrompaf.fa
```


```shell
lastz GCA_idfrompaf.fa[multiple] UCSC_idfrompaf.fa[multiple] --notransition --step=20 --nogapped --format=rdotplot --ambiguous=iupac >  algnmorechrandmorechr.txt
```
Using my [R script](dotplot.R): Rscript dotplot.R input 

![plotaln4seqvs4seq](/img/aln4seqvs4seq.png)


2. If you want to align everything against everything:
```shell
cat GCA_idfrompaf.fa UCSC_idfrompaf.fa > GCA+UCSC.fa
lastz GCA+UCSC.fa[multiple] GCA+UCSC.fa[multiple] --step=20 --nogapped  --format=maf > GCA+UCSC.maf #this output is big
```
3. If you want to extract a single region for both references

http://emboss.sourceforge.net/apps/cvs/emboss/apps/extractseq.html

The same commands for two references:

```shell
extractseq
Input sequence:human/GCA_000001405.27_GRCh38.p12_genomic_changeid.fa
Regions to extract (eg: 4-57,78-94) [1-248956422]: 30000000-60000000
output sequence(s) [chr1.fasta]:  GCAregion.fa 
```

```shell
lastz GCAregion.fa UCSCregion.fa --step=20 --nogapped --format=maf > alnsameregion.maf 
last-dotplot alnsameregion.maf sameregionshumanvsmouse.png
lastz GCAregion.fa UCSCregion.fa --step=20 --nogapped --format=rdotplot > humanvsmousesameregion.txt #forvizwithR. The output is bad, there isn't the same region syntenic share between two species
```

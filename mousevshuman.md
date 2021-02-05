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
last-dotplot aln.maf algn.png
```
or with R

```shell
lastz GCA_idfrompaf.fa UCSC_idfrompaf.fa --notransition --step=20 --nogapped --format=rdotplot > humanvsmouse.txt
```
```R
dots = read.table("aln.txt",header=T)
plot(dots,type="l")
```
or in R with my script: 
```R
library(ggplot2)
library(tidyverse)
myd <- read.table("humanvsmouse.txt", header=TRUE, sep="\t", row.names=NULL)
x <- myd[,1, drop=FALSE]
x$homo  <- "homo"
y <- myd[,2, drop=FALSE]
y$mouse  <- "mouse"
df = data.frame (x, y)
colnames(df) = c("value_homo", "id_homo", "value_mouse", "id_mouse")
myd = df %>% gather(species,values,starts_with("value_"))
p = ggplot(myd, aes(values, values)) + geom_point(aes(colour = as.factor(species)))
q = ggplot(myd, aes(species, values)) + geom_point(aes(colour = as.factor(species)))
```



I tried with another IDS: 



I tried with more IDS:




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
lastz GCA_pariamm.fa UCSC_pariamm.fa --step=20 --nogapped --format=rdotplot > humanvsmousesameregion.txt #forvizwithR
```

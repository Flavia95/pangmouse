# Mapping reference genomes (mouse against human)

**Step before starting**

- ### I changed IDs of reference genomes:

```shell
awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"Mouse#"$1} else {print $0}}' UCSC_mm10.fa | grep -o '^\S*' > UCSC_mm10.fa_changeid.fa
```
```shell
awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"Homosapiens#"$1} else {print $0}}' GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna | grep -o '^\S*' > GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_changeid.fa
```

- ### I tried lastz:

I downloaded these sequences from here (https://www.ensembl.org/Homo_sapiens/Location/Synteny?r=17:63992802-64038237;otherspecies=Mus_musculus)

PRR29 (ENSG00000224383)	17:63998351-64004305	→	Prr29 (ENSMUSG00000009210)	11:106365472-106377558

I used lastz for the alignment.

```shell
lastz chr17hum.fa chr11mouse.fa --notransition --step=20 --nogapped --format=maf > aln.maf
last-dotplot aln.maf algn.png
```
![algn.png](/img/algn.png)

## Approximate mapping

I used wfmash

```shell
wfmash -m -K -k 29 -p 75 -s 10000000 -l 0 -t 16 GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_changeid.fa.gz UCSC_mm10.fa_changeid.fa.gz > human-mouse.-m-k29-p75-K-s10000000-l0.paf
```

```shell
sort -n -k 11 human-mouse.-m-k29-p75-K-s10000000-l0.paf | awk '{ sum += $11 } END { print sum }'
1700012320
```

## Possible use of approximate mapping + lastz

1. **Taked a single line from the PAF obtained by the approximate mapping (wfmash)**

I Used the coordinates and strand orientation to extract the query and target subsequence and then I aligned these with lastz.

2. **I extracted the corresponding sequences from the fasta files**

#filename: extractseq_mousefrompaf.bed ---> Mouse#chrX	30000000	60000000	-                   

#filename: extractseq_humanfrompaf.bed --> Homosapiens#chr1	44910958	73919870	-  

```shell
bedtools getfasta -fi mouse/UCSC_mm10_changeid.fa -bed extractseq_mousefrompaf.bed > UCSC_idfrompaf.fa
bedtools getfasta -fi human/GCA_000001405.27_GRCh38.p12_genomic_changeid.fa -bed extractseq_humanfrompaf.bed > GCA_idfrompaf.fa
```

3. **I used lastz for alignment**

```shell
lastz GCA_idfrompaf.fa UCSC_idfrompaf.fa --notransition --step=20 --nogapped --format=maf > aln.maf
```
4. **I viewed the output using three different methods, the first two seems don't work very well**

- With last-dotplot:

```shell
last-dotplot aln.maf algnChrxandChr1.png
```
![algnChrxandChr1.png](/img/algnChrxandChr1.png)

- With the output obtained by lastz, using R:

```shell
lastz GCA_idfrompaf.fa UCSC_idfrompaf.fa --notransition --step=20 --nogapped --format=rdotplot >  algnChrxandChr1.txt
```
```R
dots = read.table("algnChrxandChr1.txt",header=T,row.names=NULL)
plot(dots,type="l")
```
![sample.png](/img/sample.png)

- With my script using the output obtained by Lastz (algnChrxandChr1.txt) using R:
 
[dotplot.R](script/dotplot.R) Rscript dotplot.R input 

![alnchrx_chr1.png](/img/algnChrxandChr1.png)

## I tried with more IDS:

- Using the IDs extracted by the PAF file:

[extractseq_humanfrompaf.bed](test/extractseq_humanfrompaf.bed)

[extractseqmousefrompaf.bed](test/extractseq_mousefrompaf.bed)

```shell
bedtools getfasta -fi mouse/UCSC_mm10_changeid.fa -bed extractseq_mousefrompaf.bed > UCSC_idfrompaf.fa
bedtools getfasta -fi human/GCA_000001405.27_GRCh38.p12_genomic_changeid.fa -bed extractseq_humanfrompaf.bed > GCA_idfrompaf.fa
```

- Aligment using Lastz:

```shell
lastz GCA_idfrompaf.fa[multiple] UCSC_idfrompaf.fa[multiple] --notransition --step=20 --nogapped --format=rdotplot --ambiguous=iupac >  algnmorechrandmorechr.txt
```

- Using my [R script](dotplot.R) for visualized algnmorechrandmorechr.txt

Rscript dotplot.R input 

![plotaln4seqvs4seq](/img/aln4seqvs4seq.png)

## Other considerations

1. If we want to align everything against everything:

```shell
cat GCA_idfrompaf.fa UCSC_idfrompaf.fa > GCA+UCSC.fa
lastz GCA+UCSC.fa[multiple] GCA+UCSC.fa[multiple] --step=20 --nogapped  --format=maf > GCA+UCSC.maf #this output is big
```
2. If we want to extract a single region for both references:

I used it http://emboss.sourceforge.net/apps/cvs/emboss/apps/extractseq.html

The same commands for two reference genomes:

```shell
extractseq
Input sequence:human/GCA_000001405.27_GRCh38.p12_genomic_changeid.fa
Regions to extract (eg: 4-57,78-94) [1-248956422]: 30000000-60000000
output sequence(s) [chr1.fasta]:  GCAregion.fa 
```

```shell
lastz GCAregion.fa UCSCregion.fa --step=20 --nogapped --format=maf > alnsameregion.maf 
last-dotplot alnsameregion.maf sameregionshumanvsmouse.png
lastz GCAregion.fa UCSCregion.fa --step=20 --nogapped --format=rdotplot > humanvsmousesameregion.txt #forvizwithR. 
```
I don't think this makes sense because there is no single (syntenic) region that is the same for both species, even taking into consideration a single chromosome.

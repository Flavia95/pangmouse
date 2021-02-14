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
4. **I viewed the output using different methods, but the visualization with the maf output doesn't work**

- With last-dotplot:

```shell
last-dotplot aln.maf algnChrxandChr1.png
```
![algnChrxandChr1.png](/img/algnChrxandChr1.png)


## Use approximate mapping + lastz (with paf output)

With one ID for species, seems works:

The same for more IDs:
```shell
wfmash -m -s 100000 -p 75 UCSC_mm10_chr11_only.fa GCA_chr17_only.fa > UCSCvsGCA_chr11_chr17.paf && awk -F "\t" 'OFS="\t" {print $1, $3, $4, $5 > ("extractseq_homo.bed")}' UCSCvsGCA_chr11_chr17.paf && awk -F "\t" 'OFS="\t" {print $6, $8, $9, $5 > ("extractseq_mouse.bed")}' UCSCvsGCA_chr11_chr17.paf && bedtools getfasta -fi GCA_chr17_only.fa -bed extractseq_homo.bed > GCA_idfrompafchr17_moreseq.fa && bedtools getfasta -fi UCSC_mm10_chr11_only.fa -bed extractseq_mouse.bed > UCSC_idfrompafchr11_moreseq.fa && ./lastz UCSC_idfrompafchr11_moreseq.fa[multiple] GCA_idfrompafchr17_moreseq.fa[multiple] --format=paf:wfmash --gfextend --nochain --gapped > GCA_UCSC_lastz.paf && sort -n -k 8 GCA_UCSC_lastz.paf > GCA_UCSC_lastz_sort.paf && ./paf2dotplot png large GCA_UCSC_lastz_sort.paf
```
The results are good but I used another command of lastz:


1. Convert the two genomes to 2bit (full FASTA files you used as input to wfmash)

2. For each line in the wfmash PAF

3. Run lastz using the 2bit references and the new command line:
file.2bit/chr3[100..200]

cat UCSCvsGCA_chr11_chr17.paf | while read line ; do ./lastz  --nogfextend --nochain --nogapped --format=paf:wfmash UCSC_prova.2bit/Mouse#chr11[69274866..97580505] GCA_prova.2bit/Homosapiens#chr17[7700000..38600000]; done> UCSCvsGCA_lastz_after2bit.paf






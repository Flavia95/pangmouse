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

1. Approximate mapping based on the mashmap2 distance, using wfmash:
```shell
mash sketch -s 1000000 GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna &
mash sketch -s 1000000 UCSC_mm10.fa &
wait
mash dist UCSC_mm10.fa.msh GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.msh
UCSC_mm10.fa    GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna     0.187839        0       9774/1000000
```

Setting the segment length to 10 Mbp, a very large kmer size for the mash calculation (minimum hash of kmers) and an identity threshold of at least 80% seems to give us a rough mapping that covers at least the pair of genomes.

2. All parameters tested and the total length of the mapped sequences:
```shell
human-mouse.-m-p90-s5000-l0.paf 2485356
human-mouse.-M-m-k19.p85-s200000-l0.paf 8800000
human-mouse.-M-m-k23-p85-s100000-l0.paf 287800000
human-mouse.-m-k29-p70-K-s10000000-l0.paf 1055025461
human-mouse.-M-m-k27.p85-s200000-l0.paf 1146600000
human-mouse.-m-k29-p75-K-s5000000-l0.paf 1314683097
human-mouse.-m-k29-p77-K-s10000000-l0.paf 1422166719
human-mouse.-m-k29-p77-K-s5000000-l0.paf 1437019121
human-mouse.-M-m-k29.p85-s200000-l0.paf 1447600000
human-mouse.-m-k29-p80-K-s10000000-l0.paf 1490981665
human-mouse.-m-k29-p75-K-s10000000-l0.paf 1700012320
human-mouse.-M-m-k23-p80-s200000-l0.paf 1885400000
human-mouse.-M-m-k27.p80-s200000-l0.paf 2457000000
human-mouse.-M-m-k29.p80-s200000-l0.paf 2661400000
human-mouse.-M-m-k29-p70-K-s10000000-l0.paf 2860000000
human-mouse.-M-m-k29-p75-K-s10000000-l0.paf 2860000000
human-mouse.-M-m-k29-p77-K-s10000000-l0.paf 2860000000
human-mouse.-M-m-k29-p80-K-s10000000-l0.paf 2860000000
```

3. I worked on this:

```shell
wfmash -m -K -k 29 -p 75 -s 10000000 -l 0 -t 16 GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_changeid.fa.gz UCSC_mm10.fa_changeid.fa.gz > human-mouse.-m-k29-p75-K-s10000000-l0.paf
```

![approximatemapping](/img/approximate_mappingp75.png)

```shell
sort -n -k 11 human-mouse.-m-k29-p75-K-s10000000-l0.paf | awk '{ sum += $11 } END { print sum }'  
1700012320  #number of mapped segments
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


## Use approximate mapping + lastz, fortunately, lastz now has a PAF as an output

I used chr17 of human and chr11 of mouse, because there were synthenic regions.

I reproduce the point 2 and it works:

![chr11vsChr17.png](/img/chr11vschr17.png)


The same for more IDs:
```shell
wfmash -m -s 100000 -p 75 UCSC_mm10_chr11_only.fa GCA_chr17_only.fa > UCSCvsGCA_chr11_chr17.paf && awk -F "\t" 'OFS="\t" {print $1, $3, $4, $5 > ("extractseq_homo.bed")}' UCSCvsGCA_chr11_chr17.paf && awk -F "\t" 'OFS="\t" {print $6, $8, $9, $5 > ("extractseq_mouse.bed")}' UCSCvsGCA_chr11_chr17.paf && bedtools getfasta -fi GCA_chr17_only.fa -bed extractseq_homo.bed > GCA_idfrompafchr17_moreseq.fa && bedtools getfasta -fi UCSC_mm10_chr11_only.fa -bed extractseq_mouse.bed > UCSC_idfrompafchr11_moreseq.fa && ./lastz UCSC_idfrompafchr11_moreseq.fa[multiple] GCA_idfrompafchr17_moreseq.fa[multiple] --format=paf:wfmash --gfextend --nochain --gapped > GCA_UCSC_lastz.paf && sort -n -k 8 GCA_UCSC_lastz.paf > GCA_UCSC_lastz_sort.paf && ./paf2dotplot png large GCA_UCSC_lastz_sort.paf
```
The results are good but I used another command of lastz, for:


1. Convert the two genomes to 2bit (full FASTA files you used as input to wfmash):
```shell
faToTwoBit genome.fa genome.2bit
```
2. For each line in the wfmash PAF, run lastz using the 2bit references and the new command line:
```shell
/lastz --nogfextend --nochain --nogapped --format=paf:wfmash UCSC.2bit/Mouse#chr11[69274866..69374866] GCA.2bit/Homosapiens#chr17[7700000..7800000]; done> UCSCvsGCA_lastz_after2bit.paf && sort -n -k 8 GCA_UCSC_lastz_after2bit.paf > GCA_UCSC_lastz_sort_after2bit.paf && ./paf2dotplot png large GCA_UCSC_lastz_sort_after2bit.paf
```
3. Plot the query and ref positions of each PAF record using R for all records..
![refvsquery.png](/img/ref_query.png)






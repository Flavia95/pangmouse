# Mapping  reference of mousevshuman

1. I changed ID for each sequence in each reference:

```shell
awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"Mouse#"$1} else {print $0}}' UCSC_mm10.fa | grep -o '^\S*' > UCSC_mm10.fa_changeid.fa
```
```shell
awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"Homosapiens#"$1} else {print $0}}' GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna | grep -o '^\S*' > GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_changeid.fa
```

2. wfmash

```shell
wfmash -m -K -k 29 -p 75 -s 10000000 -l 0 -t 16 GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_changeid.fa.gz UCSC_mm10.fa_changeid.fa.gz > human-mouse.-m-k29-p75-K-s10000000-l0.paf
```

```shell
sort -n -k 11 human-mouse.-m-k29-p75-K-s10000000-l0.paf | awk '{ sum += $11 } END { print sum }'
1700012320
```

```shell
bedtools getfasta -fi mouse/UCSC_mm10_changeid.fa -bed extractseq_mousefrompaf.bed > UCSC_idfrompaf.fa
bedtools getfasta -fi human/GCA_000001405.27_GRCh38.p12_genomic_changeid.fa -bed extractseq_humanfrompaf.bed > GCA_idfrompaf.fa
```

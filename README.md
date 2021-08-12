# 10x technology-pggb analysis pipeline

## Data and preprocessing

1. Apply name prefixes to the reference sequences and all samples:
```
mkdir assemblies
mkdir alignments
mkdir parts
mkdir whole
awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"REF#"$1} else {print $0}}' UCSC_mm10.fa | grep -o ' ^\S*' > UCSC_mm10_changeid.fa #the same for all samples
```

2. Partition the assembly contigs by chromosome by mapping each assembly against the scaffolded reference, and then subsetting the graph. Here we use [wfmash](https://github.com/ekg/wfmash) for the mapping.

```
cd ..
ls assemblies | grep '_supernova_changeid.fasta$' | cut -f 1-2 -d . | sort -V | uniq >haps.list
ref=assemblies/UCSC_mm10_changeid.fa
for hap in $(cat haps.list);
do
    in=assemblies/$hap
    out=alignments/$hap.vs.ref.paf
    wfmash -t 25 -m -N -p 90 -s 50000 $ref $in >$out
done

```

3. Subset by chromosome.

```
(seq 19; echo X; echo Y; echo M) | while read i; do awk '$6 ~ "chr'$i'$"' $(ls alignments/*.vs.ref.paf | sort -V) | cut -f 1 | sort -V >parts/chr$i.contigs; done
cat assemblies/*.fa > 148strains.fa
(seq 19;echo X;echo Y;echo M) | while read i; do xargs samtools faidx assemblies/148strains.fa  < parts/chr$i.contigs > parts/chr$i.pan.fa; done
```

4. Then we can merge the reference scaffold and contigs:

```
(seq 19; echo X; echo Y; echo M) | while read i; do samtools faidx UCSC_mm10_changeid.fa REF#chr$i > UCSC_mm10_changeid.chr$i.fa && cat UCSC_mm10_changeid.chr$i.fa parts/chr$i.pan.fa > chr$i.pan+ref.fa; done
```

We will use these files directly in [pggb](https://github.com/pangenome/pggb).


## Graph generation

We apply [pggb](https://github.com/pangenome/pggb) and variant calling, here an example for chr19.

```
sbatch -p lowmem -c 48 --wrap 'cd /scratch && pggb -t 48 -i /lizardfs/flaviav/mouse/chr19.pan+ref.fa.gz -s 100000 -p 98  -n 140 -k 229  -B 10000000 -w 1000000 -G 13219,15331 -T 24 -P 1,19,39,3,81,1 -v -L  -V REF:/lizardfs/flaviav/out/chr19.pan.mouse/sample.names -o /scratch/chr19.pan.mouse; mv /scratch/chr19.pan.mouse '$(pwd)
```
Adjust pangenomic VCF with [adjustpangvcf.sh](script/adjustpangvcf.sh)

## Graph evaluation
```
mkdir evaluation
```
Prepare the reference genome:
```
rtg format -o ref.sdf /home/davida/UCSC_mm10.fa
```
Prepare the truth genomic set:

[adjustvcf.sh](script/adjustvcf.sh) genomic.vcf

Run the evaluations. 
Run:https://github.com/pangenome/HPRCyear1v2genbank/blob/main/evaluation/vcf_evaluation.sh


## Pangenome of chr19 from nanopore reads of 2 strains (C57BL6J, DBA2J)
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

## 1. Correction of nanopore reads

*flaviav@penguin2:/home/flaviav/data/nanopore_parental*

- Conversion supernova assembly to graph
```shell
MBG -i /home/flaviav/data/nanopore_parental/Supernovadenovo_raw/DBA2J_supernova_3lanes_de_novo_raw.changeid.fasta -o /home/flaviav/data/nanopore_parental/MBG+GA/DBA2J.gfa -k 1667 -w 1447 -a 1 -u 1 -t 16 --no-hpc

MBG -i /home/flaviav/data/nanopore_parental/Supernovadenovo_raw/C57BL6J_supernova_2lanes_de_novo_raw.changeid.fasta -o /home/flaviav/data/nanopore_parental/MBG+GA/C57BL6J.gfa -k 1667 -w 1447 -a 1 -u 1 -t 16 --no-hpc
```

Graphs           | Nodes |Edges | Assembly size | N50       
--------------| ----|---------|---------|------
DBA2J (3.6 G) | 588,196 |  78,659 |3,767,966,689 bp | 8,484 |
C57BL6J (4,2 G) |831,236 | 387,727 | 4,383,396,388 bp | 7,094 |

- Aligned nanopore reads against graphs (obtained from the previous step)

*flaviav@penguin2:/home/flaviav/data/nanopore_parental/MBG+GA*
```shell
GraphAligner -g DBA2J.gfa -f /home/flaviav/data/nanopore_parental/D/home/flaviav/data/nanopore_parental/nanopore_raw/DBA2J.changeid.fa --corrected-out DBA2J_corrected.fa -x dbg -t 3

GraphAligner -g C57BL6J.gfa -f /home/flaviav/data/nanopore_parental/D/home/flaviav/data/nanopore_parental/nanopore_raw/C57BL_6J.changeid.fa --corrected-out C57BL6J_corrected.fa -x dbg -t 3
```
## 2. Evaluation of corrected reads

- Mapped corrected and raw nanopore reads against whole reference genome

*flaviav@penguin2:/home/flaviav/data/nanopore_parental/winnowmap*
```shell
- winnowmap -W repetitive_k15.txt -ax map-ont UCSC_mm10.fa DBA2J_corrected.fa >DBA2J_corrected_winnow.sam
- winnowmap -W repetitive_k15.txt -ax map-ont UCSC_mm10.fa C57BL6J_corrected.fa >C57BL6J_corrected_winnow.sam
- winnowmap -W repetitive_k15.txt -ax map-ont UCSC_mm10.fa nanopore_raw/DBA2J.changeid.fa > DBA2J_uncorrected_winnow.sam
- winnowmap -W repetitive_k15.txt -ax map-ont UCSC_mm10.fa nanopore_raw/C57BL_6J.changeid.fa > C57BL6J_uncorrected_winnow.sam
```
samtools stats: *flaviav@penguin2:/home/flaviav/data/nanopore_parental/stats_correctedreads*

## 3. Assembly with CANU



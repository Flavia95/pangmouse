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
DBA2J |  |  
C57BL6J | | 


GRAPH FOR DBA2J:
nodes: 588,196
edges: 78,659
assembly size 3,767,966,689 bp, N50 8484
3.6 G
GRAPH FOR C57BL6J:
nodes: 831,236
edges: 387,727
assembly size 4,383,396,388 bp, N50 7094



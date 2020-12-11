1. Change id for each strains

awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"BXD001#"$1} else {print $0}}' /home/flaviav/data/BXD001_supernova_de_novo.fasta | grep -o '^\S*' > BXD001_supernova_changeid.fasta


Variants        | Method pangenome (call variants: all vs all)  | Method pangenome (call variants: ref vs all) | Method standard
--------------| ------------- | --------- | ------                                                     
no. of SNP |        484  |   165  |  101672                                           
no. of MNP|        25     |     9  |  ----            
no. of INDEL|         69    | 25   |  19853
no. of SNP/MNP  |          4 |  2  |  ----
no. of SNP/INDEL|         76  |  26 |  ----
no. of MNP/INDEL |         69  |   23  | ----
no. of MNP/CLUMPED|         12  |    4 | ----
no. of INDEL/CLUMPED|   341 | 114 |  ----
no. of micro variants  |       1080 |   368 |  121525
no. of clumped variants |       353  |  118 | ----
no. of block substitutions |         41 |   15 | ----
no. of complex substitutions |        486 |   163  |  ----
no. of VCF records |       1080   |  368 |  121525
                                                              

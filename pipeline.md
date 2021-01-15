## Build pangenome for 10 strains (C57BL6J+DBA2J+BXD001+BXD002+BXD005+BXD006+BXD008+BXD009+BXD011+BXD012)

1. For simplicity I changed ID for each strain

awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"DBA2J#"$1} else {print $0}}' /home/flaviav/data/DBA2J_supernova_3lanes_de_novo.fasta | grep -o '^\S*' > DBA2J_supernova_changeid.fasta


       

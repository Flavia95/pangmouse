1. Change id for each strains

awk -F "::" '{if($1~">"){gsub(">","");print ">"$2"BXD001#"$1} else {print $0}}' /home/flaviav/data/BXD001_supernova_de_novo.fasta | grep -o '^\S*' > BXD001_supernova_changeid.fasta

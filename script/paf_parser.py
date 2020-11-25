import sys

strain_paf=sys.argv[1]
strain_chr19_paf= sys.argv[2]
id_strain = sys.argv[3]

with open(strain_paf, 'r') as f, open(strain_chr19_paf, 'w') as paf, open(id_strain, 'w') as id:
    for line in f:
        if 'chr19' in line:
            paf.write(line)
            x = (line.split('\t')[0])
            id.write(str(x)+ '\n')

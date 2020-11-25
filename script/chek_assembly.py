import sys

fasta = sys.argv[1]

def fastaParser(infile):
    seqs = []
    headers = []
    with open(infile) as f:
        sequence = ""
        header = None
        for line in f:
            if line.startswith('>'):
                headers.append(line[1:-1])
                if header:
                    seqs.append([sequence])
                sequence = ""
                header = line[1:]
            else:
                sequence += line.rstrip()
        seqs.append([sequence])
    return headers, seqs

headers, seqs = fastaParser(fasta)

flat_seqs = [item for sublist in seqs for item in sublist]

h = (len(headers))  #number of header = number of sequences
#print("totalnumberofsequences:",h)


def countNucs(instring):
    # will count upper and lower case sequences, if do not want lower case remove .upper()
    g = instring.upper().count('G')
    c = instring.upper().count('C')
    a = instring.upper().count('A')
    t = instring.upper().count('T')
    n = instring.upper().count('N')  #for gap
    x = (g+c+a+t+n) #lenofsequences


    return 'G = {}, C = {}, A = {}, T = {}, N = {}, lenseq:{}'.format(g, c, a, t,n,x)

lista = []

for header, seq in zip(headers, flat_seqs):
    lista.append([header, countNucs(seq)])
#print(lista)


with open('DBA2J_supernova_de_novo_statistics.txt', 'w') as fw:
    fw.write('##totatalnumberofsequences=' + str(h) + '\n')
    fw.write('\t'.join(['#SEQUENCES', 'COUNT_NUCLEOTIDE', 'LENSEQUENCE']) + '\n')
    fw.write('\n'.join(['\t'.join(row) for row in lista]) + '\n')

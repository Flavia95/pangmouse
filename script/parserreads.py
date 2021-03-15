from itertools import islice 
import re, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-input',type=str, help='path to the sam file', required=True)
parser.add_argument('-output',type=str, help ='path to output file txt ', required=True)
args = parser.parse_args()

idreadscentromer = []
with open(args.input) as f:
    for line in f:
            reads_id = line.split("\t")[0]
            position_left= line.split("\t")[3]
            seqlength = len(line.split("\t")[9])

            if int(position_left) >= seqlength/0.5:
                print(False)
            else:
                idreadscentromer.append(reads_id)

with open(args.output, 'w') as fw:
    for x in idreadscentromer:
            fw.write(str(x)+ '\n')

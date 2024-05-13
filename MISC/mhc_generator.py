import glob

''
fasta_file=glob.glob('input/*fasta')[0]
sequence=[]
with open(fasta_file, 'r') as fasta:
    for line in fasta:
        if line.startswith('>'):
            continue
        else:
            sequence=line.strip()


with open('input/input_mhc.txt', 'w') as mhc:
    for i in range(len(sequence)-9):
        mhc.write(f'{sequence[i:i+9]}\n')
        

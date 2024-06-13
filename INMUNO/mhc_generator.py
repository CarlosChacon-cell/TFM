import glob

''
fasta_files=glob.glob('fastas/*fasta')


for fasta_file in fasta_files:
    sequence=[]
    fasta_name=fasta_file.split('/')[1].split('.')[0]
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                continue
            else:
                sequence=line.strip()


    with open(f'input/{fasta_name}.txt', 'w') as mhc:
        for i in range(len(sequence)-9):
            mhc.write(f'{sequence[i:i+9]}\n')
        

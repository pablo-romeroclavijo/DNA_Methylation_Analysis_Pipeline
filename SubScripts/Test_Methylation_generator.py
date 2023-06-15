import random
import os
import shutil

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))

def random_methylation(data, filename):

    if not os.path.exists(directory):
        os.mkdir(directory)

    with open (data) as fp:
        with open('00_Random_methylated/' + filename, 'w') as mp:
            for name, seq in read_fasta(fp):
                methylation_probability = random.random()
                cp_positions = []

                new_seq = list(seq)
                for j in range(6, len(new_seq) - 1):
                    if new_seq[j] == 'C':
                        if new_seq[j + 1] == 'G':
                            cp_positions.append(j)

                new_seq = list(seq)
                for j in range(6, len(new_seq) - 1):
                    if new_seq[j] == 'C':
                        if j in cp_positions:
                            if random.random() <= methylation_probability:    ##methylation of CpG
                                new_seq[j] = 'T'
                        else:                                                 ## uncomplete conversion
                            if random.random() <= 0.1:
                                new_seq[j] = 'C'
                            else:
                                new_seq[j] = 'T'

                for j in range(len(new_seq) - 1):                           ## random errors
                        if random.random() <= 0.005:
                            bases =['A', 'G', 'T', 'C']
                            new_seq[j] = random.choice(bases)

                new_seq = ''.join(new_seq)
                mp.write(name)
                mp.write('\n')
                mp.write(new_seq)
                mp.write('\n')

os.chdir('/Users/romerop/Desktop/AAA_PhD/Experiments/Methylation/Analysis/Decoder/')
directory = '00_Random_methylated'
if os.path.exists(directory):
    shutil.rmtree(directory)


source = '000_Random_pre_methylated'
for filename in os.listdir(source):
    data = os.path.join(source, filename)
    random_methylation(data, filename)


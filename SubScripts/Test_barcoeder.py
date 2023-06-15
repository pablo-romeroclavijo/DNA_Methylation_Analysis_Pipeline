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


def barcode_dict():
    dict = {}
    with open('Barcodes.txt') as barcodes:
        for name, seq in read_fasta(barcodes):
            dict[name]= seq
    return dict

##indexing
directory = '000_Random_pre_methylated/'
if os.path.exists(directory):
    shutil.rmtree(directory)

for i in range(1, 3):
        with open('0_Reference/Ref_uncon.txt') as wp:
            barcodes = barcode_dict()
            print(barcodes)

            if not os.path.exists(directory):
                os.makedirs(directory)

            file = directory + 'Sample' + str(i)
            with open(file, 'w') as fp:
                for name, seq in read_fasta(wp):
                    for j in range (1, 1001):
                        choice_name = random.choice(list(barcodes.keys()))
                        choice_seq = barcodes[choice_name]

                        indexed_name = '>Seq' + str(j) + '__' + choice_name + ':' + choice_seq
                        indexed_seq = choice_seq + seq

                        fp.write(indexed_name)
                        fp.write('\n')
                        fp.write(indexed_seq)
                        fp.write('\n')
                fp.close()

        with open(file) as fp:
            for name, seq in read_fasta(fp):
                print(name,seq)






import os
import shutil
from Bio import SeqIO

def barcode_dict():
    dict = {}
    with open('0_Reference/Barcodes.txt') as barcodes:
        for record in SeqIO.parse(barcodes, 'fasta'):
            dict[record.id] = record.seq

    barcodes.close()
    return dict


def demultiplex(sample, filename):
    with open(sample) as fp:
        barcodes = barcode_dict()
        for record in SeqIO.parse(fp, 'fasta'):
            start = ''.join(list(record.seq)[0:6])

            if start in barcodes.values():
                for barcode in barcodes.items():
                    if start == barcode[1]:
                        record.id = record.id + barcode[0] + ':' + barcode[1]
                        record.seq = record.seq[6:]
                        file_name = str('1_Decoded_data/' + filename + '_' + barcode[0] + '.txt')

                        with open(file_name, 'a') as fp:
                            fp.write('>' + str(record.id))
                            fp.write('\n')
                            fp.write(str(record.seq))
                            fp.write('\n')
            else:
                new_name = record.id + "__" + ''.join(record.seq[0:6])
                record.seq = record.seq[6:]
                file_name_2 = '1_Decoded_data/' + filename + '_unmatched.txt'

                with open(file_name_2, 'a') as fp:
                    fp.write('>' + str(new_name))
                    fp.write('\n')
                    fp.write(str(record.seq))
                    fp.write('\n')

    fp.close()



working_directory = input("Choose working directory: ")
if not working_directory:
    working_directory = os.getcwd()
else:
    os.chdir(working_directory)
print("working_directory:", working_directory)

if not os.path.isfile(working_directory + '/0_Reference/Barcodes.txt'):
    print('No Barcode.txt file found')

if not os.path.isdir(working_directory + '/0_Raw_data'):
    print('No 0_Raw_data directory')

else:
    raw_data = '0_Raw_data'
    if os.path.exists("1_Decoded_data"):
        shutil.rmtree('1_Decoded_data')

    os.mkdir('1_Decoded_data')

    omitted_files = []
    for filename in os.listdir(raw_data):
        data = os.path.join(raw_data, filename)
        try:
            demultiplex(data, filename)
            print('Processed file: ', filename)
        except:
            omitted_files.append(filename)

    if omitted_files:
        print('Omitted files: ', omitted_files)
    else:
        print('All files processed')

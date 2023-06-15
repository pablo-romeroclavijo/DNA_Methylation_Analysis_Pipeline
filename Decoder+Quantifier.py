import os
import shutil
from Bio import SeqIO
from Bio import Align
import ast
import matplotlib.pyplot as plt


class Decoder:

    def __init__(self, working_directory):
        self.working_directory = working_directory
        os.chdir(self.working_directory)

    def barcode_dict(self):
        dict = {}
        with open('0_Reference/Barcodes.txt') as barcodes:
            for record in SeqIO.parse(barcodes, 'fasta'):
                dict[record.id] = record.seq

        barcodes.close()
        return dict

    def low_quality(self, seq):
        count = 0
        for n in seq:
            if n == "N":
                count += 1
        per_N = count/len(seq)
        if per_N >= 0.05:
            return True

    def check_target_seq(self, seq2):
        ref_file = '0_Reference/Ref_uncon.txt'
        with open(ref_file, 'r') as mp:
            for record in SeqIO.parse(mp, 'fasta'):
                seq1 = record.seq
        mp.close()

        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.gap_score = -3
        aligner.end_extend_gap_score = -1
        alignments = aligner.align(seq1, seq2)
        for line in alignments[0]:
            bases = list(line)
            for base in bases:
                if base == '-':
                    return True

    def demultiplex(self, sample, filename):
        with open(sample) as fp:
            barcodes = Decoder.barcode_dict(self)
            for record in SeqIO.parse(fp, 'fastq'):
                start = ''.join(list(record.seq)[0:6])

                if start in barcodes.values():
                    for barcode in barcodes.items():
                        if start == barcode[1]:
                            record.id = record.id + barcode[0] + ':' + barcode[1]
                            original_seq = list(record.seq)
                            trimmed_seq = ''.join(list(original_seq)[6:])

                            if self.low_quality(original_seq) == True:   ##removes low quality
                                new_name = record.id + "__" + barcode[0]
                                file_name_2 = '1_Decoded_data/' + filename + '_lowquality.txt'

                                with open(file_name_2, 'a') as fp:
                                    fp.write('>' + str(new_name))
                                    fp.write('\n')
                                    fp.write(str(trimmed_seq))
                                    fp.write('\n')

                            if self.check_target_seq(trimmed_seq) == True:   ##removes seq not aligning with the targer
                                new_name = record.id + "__" + barcode[0]
                                file_name_3 = '1_Decoded_data/' + filename + '_indel.txt'

                                with open(file_name_3, 'a') as fp:
                                    fp.write('>' + str(new_name))
                                    fp.write('\n')
                                    fp.write(str(trimmed_seq))
                                    fp.write('\n')

                            else:
                                file_name = str('1_Decoded_data/' + filename + '_' + barcode[0] + '.txt')
                                with open(file_name, 'a') as fp:
                                    fp.write('>' + str(record.id))
                                    fp.write('\n')
                                    fp.write(str(trimmed_seq))
                                    fp.write('\n')
                else:
                    new_name = record.id + "__" + ''.join(record.seq[0:6])
                    original_seq = list(record.seq)
                    trimmed_seq = ''.join(list(original_seq)[6:])
                    file_name_2 = '1_Decoded_data/' + filename + '_unmatched.txt'

                    with open(file_name_2, 'a') as fp:
                        fp.write('>' + str(new_name))
                        fp.write('\n')
                        fp.write(str(trimmed_seq))
                        fp.write('\n')

        fp.close()

    def run(self):

        print("working_directory:", self.working_directory)

        if not os.path.isfile(self.working_directory + '/0_Reference/Barcodes.txt'):
            print('No Barcode.txt file found')

        if not os.path.isdir(self.working_directory + '/0_Raw_data'):
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
                    Decoder.demultiplex(self, data, filename)
                    print('Processed file: ', filename)
                except:
                    omitted_files.append(filename)

            if omitted_files:
                print('Omitted files: ', omitted_files)
            else:
                print('All files processed')

class Quantifier:
    def __init__(self, working_directory):
        self.working_directory = working_directory
        os.chdir(self.working_directory)

    def rounding(self, value):
        result = round(value, 3)
        return result

    def cp_finder(self,seq_file):
        with open (seq_file) as fp:
            name = fp.readline()
            seq = fp.readline()
            new_seq = list(seq)

            cp_positions = []
            non_cp_position = []
            count = new_seq.count('C')

            for j in range(len(new_seq) - 1):
                if j != 92:
                    if new_seq[j] == 'C':
                        if new_seq[j + 1] == 'G':
                            cp_positions.append(j)
                        else:
                            non_cp_position.append(j)
            fp.close()
        return cp_positions, non_cp_position

    def quantum(self, seq, cp_positions, non_cp_position):
        new_seq = list(seq)

        methylated = []
        unmethylated = []
        uncomplete = []
        cp_errors = []
        non_cp_erros = []

        for cp in cp_positions:
            if new_seq[cp] == 'C':
                methylated.append(cp)
            elif new_seq[cp] == 'T':
                unmethylated.append(cp)
            else:
                cp_errors.append(cp)

        for non_cp in non_cp_position:
            if new_seq[non_cp] == 'C':
                uncomplete.append(non_cp)
            elif new_seq[non_cp] == 'T':
                pass
            else:
                non_cp_erros.append(non_cp)

        key_list = ["methylated", "unmethylated", "uncomplete", "cp_errors", "non_cp_erros"]
        value_list = [methylated, unmethylated, uncomplete, cp_errors, non_cp_erros]
        quantum_result = dict(zip(key_list,value_list))

        return quantum_result

    def calculation(self, data, filename, cp_positions, non_cp_position):
        if not os.path.exists("2_Methylation_results/CpG_raw-counts/"):         ##directory creation
            os.mkdir("2_Methylation_results/CpG_raw-counts/")

        if not os.path.exists("2_Methylation_results/CpG+Seq/"):         ##directory creation
            os.mkdir("2_Methylation_results/CpG+Seq/")

        if not os.path.exists("2_Methylation_results/final/"):  ##directory creation
            os.mkdir("2_Methylation_results/final/")

        file = '2_Methylation_results/CpG_raw-counts/' + filename
        file2 = '2_Methylation_results/CpG+Seq/' + filename
        file3 = '2_Methylation_results/final/' + filename


        with open(data) as mp:                              ## writes CpG-raw and CpG+Seq files
            for record in SeqIO.parse(mp, 'fasta'):
                name, seq = record.id, record.seq
                result = self.quantum(seq, cp_positions, non_cp_position)

                with open(file, 'a') as wp:
                    wp.write('>' + str(name))
                    wp.write('\n')
                    wp.write(str(result))
                    wp.write('\n')
                    wp.close()

                with open(file2, 'a') as wp:
                    name = name + ' | ' + str(result)
                    wp.write('>' + str(name))
                    wp.write('\n')
                    wp.write(str(seq))
                    wp.write('\n')
                    wp.close()
            mp.close()

        with open(file, 'r') as wp:                     ### counts per site
            q_meth = {key: None for key in cp_positions}
            q_unmeth = {key: None for key in cp_positions}
            errors = {"uncomplete":0, "cp_errors":0, "non_cp_erros":0}
            sequence_count = 0

            for record in SeqIO.parse(wp, 'fasta'):
                name, seq = record.id, record.seq
                sequence_count += 1
                result_dict = ast.literal_eval(str(seq))   ##literal evaluation to create dict form string

                sites = result_dict['methylated']           ##methylated evaluation
                for site in sites:
                    if q_meth[site] == None:
                        q_meth[site] = 1
                    else:
                        q_meth[site] += 1
                for value in q_meth.items():
                    if value[1] == None:
                        q_meth[value[0]] = 0

                sites_unm = result_dict['unmethylated'] ##unmethylated evaluation
                for site in sites_unm:
                    if q_unmeth[site] == None:
                        q_unmeth[site] = 1
                    else:
                        q_unmeth[site] += 1
                for value in q_unmeth.items():
                    if value[1] == None:
                        q_unmeth[value[0]] = 0

                for error_type in errors.keys():    ##error counting
                    for error in result_dict[error_type]:
                        errors[error_type] +=1




            per_Q_meth = q_unmeth.copy()            ##percentage calculation
            per_Q_meth.update((x, self.rounding(y/sequence_count)) for x,y in q_meth.items())

            per_Q_unmeth = q_unmeth.copy()          ##percentage calculation
            per_Q_unmeth.update((x, self.rounding(y/sequence_count)) for x, y in q_unmeth.items())

            per_errors = errors.copy()              ##percentage errors
            per_errors['uncomplete'] = self.rounding(per_errors['uncomplete']/(len(non_cp_position)*sequence_count))
            per_errors['cp_errors'] = self.rounding(per_errors['cp_errors']/(sequence_count*len(cp_positions)))
            per_errors['non_cp_erros'] = self.rounding(per_errors['non_cp_erros']/(sequence_count*len(non_cp_position)))

            correction_dict = per_Q_meth.copy()
            correction_dict.update(((x, self.rounding(y*per_errors['uncomplete'])) for x,y in correction_dict.items()))

            adj_per_meth = per_Q_meth.copy()
            adj_per_meth.update(((x, y / (1 - per_errors['cp_errors'])) for x, y in adj_per_meth.items()))          ##adjust by cp_error to make it 100%
            adj_per_unmeth = per_Q_unmeth.copy()
            adj_per_unmeth.update(((x, y / (1 - per_errors['cp_errors'])) for x, y in adj_per_unmeth.items()))      ##adjust by cp_error to make it 100%

            for x, y in correction_dict.items():
                adj_per_meth[x] = self.rounding(adj_per_meth[x] - y)         ##adjust by uncompelte conversion
                adj_per_unmeth[x] = self.rounding(adj_per_unmeth[x] + y)     ##adjust by uncompelte conversion



            overall_meth = self.rounding(sum(q_meth.values()) / (sequence_count*len(cp_positions)))
            adj_overall_meth = self.rounding(sum(adj_per_meth.values()) /(len(cp_positions)))

            overall_unmeth = self.rounding(sum(q_unmeth.values()) / (sequence_count*len(cp_positions)))
            adj_overall_unmeth = self.rounding(sum(adj_per_unmeth.values()) / (len(cp_positions)))

            overall_cp_errors = self.rounding(per_errors['cp_errors'])

            total = overall_cp_errors+overall_unmeth+overall_meth
            adj_total = adj_overall_meth+adj_overall_unmeth

            ##writes final file
            with open(file3, 'w') as fp:                ##writes final file
                fp.write(filename)
                fp.write('\n')
                fp.write("sequence_count: ")
                fp.write(str(sequence_count))
                fp.write('\n')
                fp.write("q_meth: ")
                fp.write(str(q_meth))
                fp.write('\n')
                fp.write("q_unmeth: ")
                fp.write(str(q_unmeth))
                fp.write('\n')
                fp.write('errors: ')
                fp.write(str(errors))
                fp.write('\n')
                fp.write('--------------------------------')
                fp.write('\n')
                fp.write('per_Q_meth: ')
                fp.write(str(per_Q_meth))
                fp.write('\n')
                fp.write('\n')
                fp.write('per_Q_unmeth: ')
                fp.write(str(per_Q_unmeth))
                fp.write('\n')
                fp.write('\n')
                fp.write('per_errors: ')
                fp.write(str(per_errors))
                fp.write('\n')
                fp.write('--------------------------------')
                fp.write('\n')
                fp.write("overall_meth: ")
                fp.write(str(overall_meth))
                fp.write('\n')
                fp.write("overall_unmeth: ")
                fp.write(str(overall_unmeth))
                fp.write('\n')
                fp.write('overall_errores: ')
                fp.write(str(overall_cp_errors))
                fp.write('\n')
                fp.write(str(total))
                fp.write('\n')
                fp.write('--------------------------------')
                fp.write('\n')
                fp.write('adjusted_per_Q_meth: ')
                fp.write(str(adj_per_meth))
                fp.write('\n')
                fp.write('\n')
                fp.write('adjusted_per_Q_unmeth: ')
                fp.write(str(adj_per_unmeth))
                fp.write('\n')
                fp.write('\n')
                fp.write("adjusted_overall_meth: ")
                fp.write(str(adj_overall_meth))
                fp.write('\n')
                fp.write('\n')
                fp.write("adjusted_overall_unmeth: ")
                fp.write(str(adj_overall_unmeth))
                fp.write('\n')
                fp.write('\n')
                fp.write(str(adj_total))
                fp.write('\n')
                fp.write('--------------------------------')
                fp.write('--------------------------------')
                fp.write('--------------------------------')
                fp.write('\n')
                fp.write('\n')     #wwrr


            print(filename)
            print("sequence_count: ", sequence_count)
            print("q_meth: ", q_meth)
            print("q_unmeth: ", q_unmeth)
            print('errors: ', errors)
            print('--------------------------------')
            print('per_Q_meth: ', per_Q_meth)
            print('per_Q_unmeth: ', per_Q_unmeth)
            print('per_errors: ', per_errors)
            print('--------------------------------')
            print("overall_meth: ", overall_meth)
            print("overall_unmeth: ", overall_unmeth)
            print('overall_errores: ',overall_cp_errors)
            print(total)
            print('--------------------------------')
            print('adjusted_per_Q_meth: ', adj_per_meth)
            print('adjusted_per_Q_unmeth: ', adj_per_unmeth)
            print("adjusted_overall_meth: ", adj_overall_meth)
            print("adjusted_overall_unmeth: ", adj_overall_unmeth)
            print(adj_total)
            print('\n')
            print('--------------------------------')
            print('--------------------------------')
            print('--------------------------------')
            print('\n')
            print('\n')

    def run(self):
        print("working_directory:", self.working_directory)

        ref_seq = '0_Reference/Ref_uncon.txt'               ###change reference sequence
        print("Reference_sequence: ", ref_seq, '\n -----------------\n')
        directory = '1_Decoded_data'                        ###change for source data directory

        if  os.path.exists("2_Methylation_results"):
            shutil.rmtree("2_Methylation_results")

        os.mkdir("2_Methylation_results")

        if  os.path.exists("3_Plots"):
            shutil.rmtree("3_Plots")

        os.mkdir("3_Plots")

        with open(ref_seq) as mp:                           ##create site map
            cp_positions, non_cp_position = Quantifier.cp_finder(self, ref_seq)
            with open("2_Methylation_results/AA_Reference_sequence_Analysis", 'w') as fp:
                fp.write('CpG_sites: ')
                fp.write(str(cp_positions))
                fp.write('\n')
                fp.write('\n')
                fp.write('non_CpG_sites: ')
                fp.write(str(non_cp_position))

            fp.close()
            mp.close()

        omitted_files = []
        for filename in os.listdir(directory):
            data = os.path.join(directory, filename)
            try:
                self.calculation(data, filename, cp_positions, non_cp_position)
            except:
                omitted_files.append(filename)
        if omitted_files:
            print('Omitted files: ', omitted_files)
        else:
            print('All files processed')

class Plotting:
    def  __init__(self, working_directory, queries, qc_queries, barcodes):
        self.working_directory = working_directory
        self.queries = queries
        self.qc_queries = qc_queries
        self.barcodes = barcodes
        os.chdir(self.working_directory)

    def find_data(self, query, file):
            for line in file.readlines():
                if line.startswith(query):
                    return line

    def plot(self, dictionary, query, max, y_label, filename, sample):
        if not os.path.exists("3_Plots/" + query[0]):         ##directory creation
            os.mkdir("3_Plots/" + query[0])

        file = "3_Plots/" + query[0] + '/' + filename

        x_values = []
        y_values = []
        for key, value in dictionary.items():
            x_values.append(str(key))
            y_values.append(value)

        plt.bar(x_values, y_values)
        plt.xticks(label=x_values, rotation = 90)
        plt.xlabel('CpG postision')

        plt.ylim(0, max)
        plt.ylabel(y_label)
        plt.title(query[0] +'__' + sample)

        plt.savefig(file + '.pdf')
        plt.show()

    def plot_QC(self, query, x_values, y_values):

        file = "3_Plots/" + query[1]

        plt.bar(x_values, [float(y) for y in y_values])
        plt.xticks(label=x_values, rotation=90)
        plt.xlabel('Sample')

        maximum =  query[1]
        if maximum == True:
            plt.ylim(0, maximum)
            print('maximum applied')
        plt.ylabel(query[3])
        plt.title(query[1])

        plt.savefig(file + '.pdf')
        plt.show()

    def run(self):
    ##finding data for queries
        for filename in os.listdir('2_Methylation_results/final'):
            file = '2_Methylation_results/final/' + filename

            for barcode in self.barcodes:
                sample = 'Barcode' + str(barcode)

                if sample + '.txt' in filename:

                    for query in self.queries:
                        with open(file, 'r') as fp:
                            result = Plotting.find_data(self, query[0], fp)
                            trimmed_result = result.replace(query[0] + ': ', '')
                            final_result = ast.literal_eval(trimmed_result)

                            max = query [1]
                            y_label = query[2]
                            self.plot(final_result, query, max, y_label, filename, sample)
                        fp.close()

    ##finding data for qc_queries
        for query in self.qc_queries:
            x_values = []
            y_values = []
            barcodes = query[4]
            for barcode in barcodes:
                sample = str(barcode)
                for filename in os.listdir('2_Methylation_results/final'):
                    if sample + '.txt' in filename:
                        file = '2_Methylation_results/final/' + filename

                        with open(file, 'r') as fp:
                            result = Plotting.find_data(self, query[0], fp)
                            trimmed_result = result.replace(query[0] + ': ', '')
                            trimmed_result_2 = trimmed_result.replace('\n', '')


                            x_values.append(sample)
                            y_values.append(trimmed_result_2)

                        fp.close()
            self.plot_QC(query, x_values, y_values)



##Choose working directory
working_directory = input("Choose SAMPLE directory (DEFAULT = Parent): ")
if not working_directory:
    path = os.getcwd()
    working_directory = os.path.abspath(os.path.join(path, working_directory))



## RUN Decoder
print(' ---------------------------  \n DECODER ')
decode = Decoder(working_directory)
decode.run()
print('\n ---------------------------\n')

## Run Quantification and plotting
print('\n ---------------------------  \n QUANTIFIER')
quantify = Quantifier(working_directory)
quantify.run()


print('\n ---------------------------  \n PLOTER')
queries = [('adjusted_per_Q_meth', 1, 'Methylation prevalence'), ('per_Q_meth', 1, 'Methylation prevalence'),  ('q_meth', 4000, 'Counts')]   ## queries as tuples of query, max and y_label
barcodes = [x for x in range(1,13)]

qc_queries = [('overall_meth', 'Overall_meth', 1 ,'Methylation prevalence', [('Barcode' + str(x)) for x in range(1,13)]),
              ('sequence_count', 'Sequence_counts', False, 'Counts', [('Barcode' + str(x)) for x in range(1,13)]),
              ('sequence_count', 'QC_Sequence_count', False, 'Counts', ['indel', 'lowquality', 'unmatched'])]


plotter = Plotting(working_directory, queries, qc_queries, barcodes)
plotter.run()


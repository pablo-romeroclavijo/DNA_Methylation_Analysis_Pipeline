import ast
import os
import shutil
import matplotlib.pyplot as plt
from Bio import SeqIO

def rounding(value):
    result = round(value, 3)
    return result

def cp_finder(seq_file):
    with open (seq_file) as fp:
        name = fp.readline()
        seq = fp.readline()
        new_seq = list(seq)

        cp_positions = []
        non_cp_position = []
        count = new_seq.count('C')

        for j in range(len(new_seq) - 1):
            if new_seq[j] == 'C':
                if new_seq[j + 1] == 'G':
                    cp_positions.append(j)
                else:
                    non_cp_position.append(j)
        fp.close()
    return cp_positions, non_cp_position

def quantum(seq, cp_positions, non_cp_position):
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

def plotting(dictionary, filename, max, y_label, title):
    if not os.path.exists("3_Plots/" + title):         ##directory creation
        os.mkdir("3_Plots/" + title)

    file = "3_Plots/" + title + '/' + filename

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
    plt.title(title)

    plt.savefig(file + '.png')
    plt.show()

def calculation(data, filename, cp_positions, non_cp_position):
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
            result = quantum(seq, cp_positions, non_cp_position)

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

    with open(file, 'r') as wp:   ### counts per site
        q_meth = {key: None for key in cp_positions}
        q_unmeth = {key: None for key in cp_positions}
        errors = {"uncomplete":0, "cp_errors":0, "non_cp_erros":0}
        sequence_count = 0

        for record in SeqIO.parse(wp, 'fasta'):
            name, seq = record.id, record.seq
            sequence_count += 1
            result_dict = ast.literal_eval(str(seq))   ##literal evaluation to create dict form string

            sites = result_dict['methylated']  ##methylated evaluation
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
        per_Q_meth.update((x, rounding(y/sequence_count)) for x,y in q_meth.items())

        per_Q_unmeth = q_unmeth.copy()          ##percentage calculation
        per_Q_unmeth.update((x, rounding(y/sequence_count)) for x, y in q_unmeth.items())

        per_errors = errors.copy()              ##percentage errors
        per_errors['uncomplete'] = rounding(per_errors['uncomplete']/(len(non_cp_position)*sequence_count))
        per_errors['cp_errors'] = rounding(per_errors['cp_errors']/(sequence_count*len(cp_positions)))
        per_errors['non_cp_erros'] = rounding(per_errors['non_cp_erros']/(sequence_count*len(non_cp_position)))

        correction_dict = per_Q_meth.copy()
        correction_dict.update(((x, rounding(y*per_errors['uncomplete'])) for x,y in correction_dict.items()))

        adj_per_meth = per_Q_meth.copy()
        adj_per_meth.update(((x, y / (1 - per_errors['cp_errors'])) for x, y in adj_per_meth.items()))          ##adjust by cp_error to make it 100%
        adj_per_unmeth = per_Q_unmeth.copy()
        adj_per_unmeth.update(((x, y / (1 - per_errors['cp_errors'])) for x, y in adj_per_unmeth.items()))      ##adjust by cp_error to make it 100%

        for x, y in correction_dict.items():
            adj_per_meth[x] = rounding(adj_per_meth[x] - y)         ##adjust by uncompelte conversion
            adj_per_unmeth[x] = rounding(adj_per_unmeth[x] + y)     ##adjust by uncompelte conversion



        overall_meth = rounding(sum(q_meth.values()) / (sequence_count*len(cp_positions)))
        adj_overall_meth = rounding(sum(adj_per_meth.values()) /(len(cp_positions)))

        overall_unmeth = rounding(sum(q_unmeth.values()) / (sequence_count*len(cp_positions)))
        adj_overall_unmeth = rounding(sum(adj_per_unmeth.values()) / (len(cp_positions)))

        overall_cp_errors = rounding(per_errors['cp_errors'])

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
            fp.write('per_Q_unmeth: ')
            fp.write(str(per_Q_unmeth))
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
            fp.write('adjusted_per_Q_unmeth: ')
            fp.write(str(adj_per_unmeth))
            fp.write('\n')
            fp.write("adjusted_overall_meth: ")
            fp.write(str(adj_overall_meth))
            fp.write('\n')
            fp.write("adjusted_overall_unmeth: ")
            fp.write(str(adj_overall_unmeth))
            fp.write('\n')
            fp.write(str(adj_total))
            fp.write('\n')
            fp.write('--------------------------------')
            fp.write('--------------------------------')
            fp.write('--------------------------------')
            fp.write('\n')
            fp.write('\n')     #wwrr

        ## plotting
        plotting(q_meth, filename, sequence_count, "Counts", 'Raw_Counts')
        plotting(per_Q_meth, filename, 1, "Methylation prevalence", 'Percentage')
        plotting(adj_per_meth, filename, 1, "Adjusted Methylation prevalence", 'Adjusted_Percentage')


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


ref_seq = '0_Reference/Ref_uncon.txt'  ###change reference sequence
directory = '1_Decoded_data'  ###change for source data directory

if  os.path.exists("2_Methylation_results"):
    shutil.rmtree("2_Methylation_results")

os.mkdir("2_Methylation_results")

if  os.path.exists("3_Plots"):
    shutil.rmtree("3_Plots")

os.mkdir("3_Plots")

with open(ref_seq) as mp:  ##create site map
    cp_positions, non_cp_position = cp_finder(ref_seq)
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
        calculation(data, filename, cp_positions, non_cp_position)
    except:
        omitted_files.append(filename)
if omitted_files:
    print('Omitted files: ', omitted_files)
else:
    print('All files processed')

















"""find new motifs for sequences in provided pre-miRNA."""

from weblogolib import *
import subprocess
import xml.etree.ElementTree
import numpy as np
import os


def ensure_dir(file_path):
    directory = os.path.dirname(file_path+'/')
    print(directory)
    if not os.path.exists(directory):
        os.makedirs(directory)


def process_input(data, database):
    list_of_seq = data.split()
    for seq_nr in range(0, len(list_of_seq)):
        this_ok = False
        if type(list_of_seq[seq_nr]) == bytes:
            list_of_seq[seq_nr] = list_of_seq[seq_nr].decode("utf-8")
        if list_of_seq[seq_nr][0:3] in ['hsa', 'mmu', 'ath']:
            filename = None
            if database.split('_')[0] == 'loops':
                filename = 'miR_loops_'+list_of_seq[seq_nr][0:3]
            if database.split('_')[0] == 'sh':
                filename = 'miR_sh_'+list_of_seq[seq_nr][0:3]
            if database.split('_')[0] == 'linking':
                filename = 'miR_linking_'+list_of_seq[seq_nr][0:3]
            if database.split('_')[0] == 'mature':
                filename = 'miR_mature_'+list_of_seq[seq_nr][0:3]
            with open('./mirnamotif/'+filename) as db:
                for line in db:
                    if list_of_seq[seq_nr].lower() == line.split(' ')[0]:
                        list_of_seq[seq_nr] = line.split(' ')[1][:-1]
                        this_ok = True
                        db.close()
                        break
            if this_ok is False:
                return ['Not all sequences or ids are correct\n' +
                        'Please check ' + list_of_seq[seq_nr]]
        else:
            seq = list_of_seq[seq_nr].upper()
            for i in seq:
                if i not in ['A', 'T', 'C', 'G', 'U']:
                    return ['Not all sequences or ids are correct\n' +
                            'Please check ' + str(seq)]
                else:
                    seq.replace('T', 'U')
            list_of_seq[seq_nr] = seq

    return list_of_seq


def motiffind(list_of_seq=None, direction='f', folder='', negative_file=''):
    """
    Find potential motif in provided sequences.

    input:
    list_of_seq - list of sequences or ids
    direction - only forward or also reverse
    output:
    file with proposed motifs
    """
    if list_of_seq is None:
        list_of_seq = []
    status = False
    if direction == 'f':
        direction = ''
    else:
        direction = ' -norc'
    fasta_from_list = ''
    counter = 1
    for seq in list_of_seq:
        fasta_from_list = fasta_from_list+'>seq'+str(counter)+'\n'+seq+'\n'
        counter = counter + 1
    ensure_dir(folder)
    with open(folder+'/input.fa', 'w') as input_file:
        input_file.write(fasta_from_list)
    command = '../meme_4.11.4/scripts/dreme-py3 -oc ' + \
              folder + ' -p ' + folder + '/input.fa -rna' + direction + \
              ' -n ' + negative_file + ' -e 0.1 -g 1000 -mink 2'
    subprocess.call(command.split())
    result_raw = xml.etree.ElementTree.parse(folder+'/dreme.xml').getroot()
    motifs_results = []
    for motifs in result_raw.findall('motifs'):
        for motif in motifs.findall('motif'):
            pfm = []
            for position in motif.findall('pos'):
                pfm.append([float(position.get('A')),
                            float(position.get('C')),
                            float(position.get('G')),
                            float(position.get('U'))])
            motifs_results.append({
                'seq': motif.get('seq'),
                'id': motif.get('id'),
                'len': motif.get('length'),
                'pvalue': motif.get('pvalue'),
                'logo_table': pfm})

    print(motifs_results)
    for motif in motifs_results:
        options = LogoOptions()
        options.title = motif['seq']
        options.color_scheme = colorscheme\
            .ColorScheme([colorscheme.SymbolColor("G", "orange"),
                          colorscheme.SymbolColor("TU", "red"),
                          colorscheme.SymbolColor("C", "blue"),
                          colorscheme.SymbolColor("A", "green")
                          ])
        data = LogoData.from_counts(alphabet='ACGU',
                                    counts=np.array(motif['logo_table']))
        format_logo = LogoFormat(data, options)
        fout = open(folder + '/' + str(motif['seq']) + '.eps', 'w')
        fout.write(eps_formatter(data, format_logo).decode("utf-8"))
        fout.close()
    file = open(folder + '/' + 'results.txt', 'w+')
    file.write('\nResults from mirnamotif - find motif\n\n')
    list_of_motifs = []
    for motif in motifs_results:
        file.write('\nmotif id = ' + motif['id'] + '\n')
        file.write('motif sequence = ' + motif['seq'] + '\n')
        list_of_motifs.append(folder + '/' + motif['seq'] + '.eps')
        file.write('length of the motif = ' + motif['len'] + '\n')
        file.write('pvalue = ' + motif['pvalue'] + '\n')
    file.close()
    if list_of_motifs:
        status = True
    return folder+'/results.txt', list_of_motifs, status

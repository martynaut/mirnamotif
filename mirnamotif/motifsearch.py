"""motif search for search of simple motifs in pre-miRNA."""

import re
import time
import os


def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)


def search(seq1='', seq2='', seq_type='loops', direction='f',
           overlap=True, database='hsa', email=''):
    """
    search function searches for motifs(seq1, seq2) in pre-miRNAs.

    input: seq1 - first or only motif
    seq2 - second motif
    seq_type - where to look for motif (loop or pre-miRNA)
    direction - only forward or also reverse
    overlap - if two motifs may overlap
    database - hsa, mmu, ath
    email - where to send results
    output: file with sequences
    """
    status = False

    seq1 = seq1.upper()
    seq2 = seq2.upper()

    iupac_dict = {
        'R': '[GA]{1}',
        'Y': '[UC]{1}',
        'M': '[AC]{1}',
        'K': '[GU]{1}',
        'S': '[GC]{1}',
        'W': '[AU]{1}',
        'H': '[ACU]{1}',
        'B': '[GUC]{1}',
        'V': '[GCA]{1}',
        'D': '[GUA]{1}',
        'N': '[GCUA]{1}',
    }
    seq1_reg = seq1
    seq2_reg = seq2
    seq1_reg_rev = seq1[::-1]
    seq2_reg_rev = seq2[::-1]
    seq1_reg_rtf = r'(\\cf1|\\cf2|2fc\\|1fc\\){0,4}'.join(list(seq1))
    seq1_reg_rev_rtf = r'(\\cf1|\\cf2|2fc\\|1fc\\){0,4}'.join(list(seq1[::-1]))
    seq2_reg_rtf = r'(\\cf1|\\cf2|2fc\\|1fc\\){0,4}'.join(list(seq2))
    seq2_reg_rev_rtf = r'(\\cf1|\\cf2|2fc\\|1fc\\){0,4}'.join(list(seq2[::-1]))
    for letter, reg in iupac_dict.items():
        seq1_reg = seq1_reg.replace(letter, reg)
        seq2_reg = seq2_reg.replace(letter, reg)
        seq1_reg_rev = seq1_reg_rev.replace(letter, reg)
        seq2_reg_rev = seq2_reg_rev.replace(letter, reg)
        seq1_reg_rtf = seq1_reg_rtf.replace(letter, reg)
        seq1_reg_rev_rtf = seq1_reg_rev_rtf.replace(letter, reg)
        seq2_reg_rtf = seq2_reg_rtf.replace(letter, reg)
        seq2_reg_rev_rtf = seq2_reg_rev_rtf.replace(letter, reg)

    hour = str(time.strftime("%H_%M_%S"))
    hourf = str(time.strftime("%H:%M:%S"))

    ensure_dir('./results/')

    filename = ('./results/'+email + '_' + str(time.strftime("%d_%m_%Y")) +
                '_' + hour + '.txt')
    filename2 = ('./results/'+email + '_' + str(time.strftime("%d_%m_%Y")) +
                 '_' + hour + '.rtf')
    file = open(filename, 'a')
    file2 = open(filename2, 'a')

    file2.write('{\\rtf1\\ansi\deff0{\colortbl;\\red0\green0\\blue0;\\red255\green0\\blue0;} ')

    file.write('\n' + (time.strftime("%d/%m/%Y")) + ' ' +
               hourf + '\n')
    file.write('\nfirst motif:\n')
    file.write(seq1 + '\n')
    if seq2 != '':
        file.write('second motif:\n')
        file.write(seq2 + '\n')
        if overlap is True:
            file.write('overlap: Yes\n')
        else:
            file.write('overlap: No\n')
    if seq_type == 'loops':
        file.write('loops sequence \n')
    else:
        file.write('pre-miRNA sequence\n')
    if direction == 'f':
        file.write('direction: forward\n')
    else:
        file.write('direction: forward and reverse\n')
    if database == 'hsa':
        file.write('database: Homo sapiens\n')
    elif database == 'mmu':
        file.write('database: Mus musculus\n')
    else:
        file.write('database: Arabidopsis thaliana\n')
    file.write(email + '\n' * 3)

    file2.write('\line ' + (time.strftime("%d/%m/%Y")) + ' ' +
                hourf + '\line ')
    file2.write('first motif:\line ')
    file2.write(seq1 + '\line ')
    if seq2 != '':
        file2.write('second motif:\line ')
        file2.write(seq2 + '\line ')
        if overlap is True:
            file2.write('overlap: Yes\line ')
        else:
            file2.write('overlap: No\line ')
    if seq_type == 'loops':
        file2.write('loops sequence\line ')
    else:
        file2.write('pre-miRNA sequence\line ')
    if direction == 'f':
        file2.write('direction: forward\line ')
    else:
        file2.write('direction: forward and reverse\line ')
    if database == 'hsa':
        file2.write('database: Homo sapiens\line ')
    elif database == 'mmu':
        file2.write('database: Mus musculus\line ')
    else:
        file2.write('database: Arabidopsis thaliana\line ')
    file2.write(email + '\line ' * 3)
    if seq_type in ['sh', 'loops']:
        mirnas = open('./mirnamotif/miR_' + seq_type + '_' + database)
    else:
        mirnas = open('./mirnamotif/miR_mature_' + database)
    if seq1[::-1] == seq1 and seq2[::-1] == seq2:
        direction = 'f'
        file.write('sequence(s) is/are palindromic\n')
        file2.write('sequence(s) is/are palindromic\line ')
    if seq2 == '' and direction == 'f':
        for line in mirnas:
            name = line.split(' ')[0]
            sequence = line.split(' ')[1][:-1]
            if bool(re.search(seq1_reg, sequence)):
                file.write(name + ' ' + sequence + '\n')
                file2.write(name + ' ' + re.sub(r'('+seq1_reg+r')', r'\cf2\1\cf1', sequence) + '\line ')
                status = True
    elif seq2 == '' and direction == 'fr':
        for line in mirnas:
            name = line.split(' ')[0]
            sequence = line.split(' ')[1][:-1]
            if bool(re.search(seq1_reg, sequence)) or bool(re.search(seq1_reg, sequence[::-1])):
                file.write(name + ' ' + sequence + '\n')
                file2.write(name + ' ' + re.sub(r'('+seq1_reg_rtf+r')',
                                                r'\cf2\1\cf1', re.sub(r'('+seq1_reg_rev_rtf+r')',
                                                                      r'\cf2\1\cf1', sequence)) + '\line ')
                status = True
    elif direction == 'f' and overlap is True:
        first = []
        second = []
        together = []
        for line in mirnas:
            name = line.split(' ')[0]
            sequence = line.split(' ')[1][:-1]
            if bool(re.search(seq1_reg, sequence)) and bool(re.search(seq2_reg, sequence)):
                together.append(name + ' ' + sequence + '\n')
                status = True
            elif bool(re.search(seq1_reg, sequence)):
                first.append(name + ' ' + sequence + '\n')
                status = True
            elif bool(re.search(seq2_reg, sequence)):
                second.append(name + ' ' + sequence + '\n')
                status = True
        file.write('\n\nfirst: ' + seq1 + '\n\n')
        file2.write('\line \line first: ' + seq1 + '\line \line ')
        for i in first:
            file.write(i)
            file2.write(re.sub(r'('+seq1_reg+r')', r'\cf2\1\cf1', i).replace('\n', '\line '))
        file.write('\n\nsecond: ' + seq2 + '\n\n')
        file2.write('\line \line second: ' + seq2 + '\line \line ')
        for i in second:
            file.write(i)
            file2.write(re.sub(r'('+seq2_reg+r')', r'\cf2\1\cf1', i).replace('\n', '\line '))
        file.write('\n\ntogether: \n\n')
        file2.write('\line \line together: \line \line ')
        for i in together:
            file.write(i)
            file2.write(re.sub(r'('+seq1_reg_rtf+r')', r'\cf2\1\cf1',
                               (re.sub(r'('+seq2_reg_rtf+r')', r'\cf2\1\cf1',
                                       i))).replace('\n', '\line '))
    elif direction == 'f' and overlap is False:
        first = []
        second = []
        together = []
        for line in mirnas:
            name = line.split(' ')[0]
            sequence = line.split(' ')[1][:-1]
            if bool(re.search(seq1_reg, sequence) and bool(re.search(seq2_reg, sequence))):
                f = re.compile(seq1_reg)
                s = re.compile(seq2_reg)
                where_f = f.search(sequence)
                where_s = s.search(sequence)
                if (where_f.end() <= where_s.start() or
                        where_s.end() <= where_f.start()):
                    together.append(name + ' ' + sequence + '\n')
                    status = True
            elif bool(re.search(seq1_reg, sequence)):
                first.append(name + ' ' + sequence + '\n')
                status = True
            elif bool(re.search(seq2_reg, sequence)):
                second.append(name + ' ' + sequence + '\n')
                status = True
        file.write('\n\nfirst: ' + seq1 + '\n\n')
        file2.write('\line \line first: ' + seq1 + '\line \line ')
        for i in first:
            file.write(i)
            file2.write(re.sub(r'('+seq1_reg+r')', r'\cf2\1\cf1', i).replace('\n', '\line '))
        file.write('\n\nsecond: ' + seq2 + '\n\n')
        file2.write('\line \line second: ' + seq2 + '\line \line ')
        for i in second:
            file.write(i)
            file2.write(re.sub(r'('+seq2_reg+r')', r'\cf2\1\cf1', i).replace('\n', '\line '))
        file.write('\n\ntogether: \n\n')
        file2.write('\line \line together: \line \line ')
        for i in together:
            file.write(i)
            file2.write(re.sub(r'('+seq1_reg_rtf+r')', r'\cf2\1\cf1',
                               (re.sub(r'('+seq2_reg_rtf+r')',
                                       r'\cf2\1\cf1', i))).replace('\n', '\line '))
    elif direction == 'fr' and overlap is True:
        if seq1[::-1] == seq1:
            file.write('sequence 1 is palindromic\n')
            file2.write('sequence 1 is palindromic\line ')
        if seq2[::-1] == seq2:
            file.write('sequence 2 is palindromic\n')
            file2.write('sequence 2 is palindromic\line ')
        first = []
        second = []
        together = []
        for line in mirnas:
            name = line.split(' ')[0]
            sequence = line.split(' ')[1][:-1]
            if ((bool(re.search(seq1_reg, sequence)) and bool(re.search(seq2_reg, sequence))) or
                    (bool(re.search(seq1_reg, sequence[::-1])) and bool(re.search(seq2_reg, sequence))) or
                    (bool(re.search(seq1_reg, sequence)) and bool(re.search(seq2_reg, sequence[::-1]))) or
                    (bool(re.search(seq1_reg, sequence[::-1])) and bool(re.search(seq2_reg, sequence[::-1])))):
                together.append(name + ' ' + sequence + '\n')
                status = True
            elif bool(re.search(seq1_reg, sequence)) or bool(re.search(seq1_reg, sequence[::-1])):
                first.append(name + ' ' + sequence + '\n')
                status = True
            elif bool(re.search(seq2_reg, sequence)) or bool(re.search(seq2_reg, sequence[::-1])):
                second.append(name + ' ' + sequence + '\n')
                status = True
        file.write('\n\nfirst: ' + seq1 + '\n\n')
        file2.write('\line \line first: ' + seq1 + '\line \line ')
        for i in first:
            file.write(i)
            file2.write(re.sub(r'('+seq1_reg_rtf+r')', r'\cf2\1\cf1',
                               re.sub(r'('+seq1_reg_rev_rtf+r')',
                                      r'\cf2\1\cf1', i)).replace('\n', '\line '))
        file.write('\n\nsecond: ' + seq2 + '\n\n')
        file2.write('\line \line second: ' + seq2 + '\line \line ')
        for i in second:
            file.write(i)
            file2.write(re.sub(r'('+seq2_reg_rtf+r')',
                               r'\cf2\1\cf1', re.sub(r'('+seq2_reg_rev_rtf+r')',
                                                     r'\cf2\1\cf1', i)).replace('\n', '\line '))
        file.write('\n\ntogether: \n\n')
        file2.write('\line \line together: \line \line ')
        for i in together:
            file.write(i)
            tog_first = re.sub(r'('+seq1_reg_rev_rtf+r')', r'\cf2\1\cf1', i)
            tog_second = (re.sub(r'('+seq2_reg_rev_rtf+r')', r'\cf2\1\cf1', tog_first))
            tog_third = (re.sub(r'('+seq2_reg_rtf+r')', r'\cf2\1\cf1', tog_second))
            file2.write(re.sub(r'('+seq1_reg_rtf+r')', r'\cf2\1\cf1', tog_third).replace('\n', '\line '))
    else:
        if seq1[::-1] == seq1:
            file.write('sequence 1 is palindromic\n')
            file2.write('sequence 1 is palindromic\line ')
        if seq2[::-1] == seq2:
            file.write('sequence 2 is palindromic\n')
            file2.write('sequence 2 is palindromic\line ')
        first = []
        second = []
        together = []
        for line in mirnas:
            name = line.split(' ')[0]
            sequence = line.split(' ')[1][:-1]
            if bool(re.search(seq1_reg, sequence)) and bool(re.search(seq2_reg, sequence)):
                f = re.compile(seq1_reg)
                s = re.compile(seq2_reg)
                where_f = f.search(sequence)
                where_s = s.search(sequence)
                if (where_f.end() <= where_s.start() or
                        where_s.end() <= where_f.start()):
                    together.append(name + ' ' + sequence + '\n')
                    status = True
            elif bool(re.search(seq1_reg_rev, sequence)) and bool(re.search(seq2_reg, sequence)):

                f = re.compile(seq1_reg_rev)
                s = re.compile(seq2_reg)
                where_f = f.search(sequence)
                where_s = s.search(sequence)
                if (where_f.end() <= where_s.start() or
                        where_s.end() <= where_f.start()):
                    together.append(name + ' ' + sequence + '\n')
                    status = True
            elif bool(re.search(seq1_reg, sequence)) and bool(re.search(seq2_reg_rev, sequence)):
                f = re.compile(seq1_reg)
                s = re.compile(seq2_reg_rev)
                where_f = f.search(sequence)
                where_s = s.search(sequence)
                if (where_f.end() <= where_s.start() or
                        where_s.end() <= where_f.start()):
                    together.append(name + ' ' + sequence + '\n')
                    status = True
            elif bool(re.search(seq1_reg_rev, sequence)) and bool(re.search(seq2_reg_rev, sequence)):
                f = re.compile(seq1_reg_rev)
                s = re.compile(seq2_reg_rev)
                where_f = f.search(sequence)
                where_s = s.search(sequence)
                if (where_f.end() <= where_s.start() or
                        where_s.end() <= where_f.start()):
                    together.append(name + ' ' + sequence + '\n')
                    status = True
            elif bool(re.search(seq1_reg, sequence)) or bool(re.search(seq1_reg_rev, sequence)):
                first.append(name + ' ' + sequence + '\n')
                status = True
            elif bool(re.search(seq2_reg, sequence)) or bool(re.search(seq2_reg_rev, sequence)):
                second.append(name + ' ' + sequence + '\n')
                status = True
        file.write('\n\nfirst: ' + seq1 + '\n\n')
        file2.write('\line \line first: ' + seq1 + '\line \line ')
        for i in first:
            file.write(i)
            file2.write(re.sub(r'('+seq1_reg_rtf+r')',
                               r'\cf2\1\cf1', re.sub(r'('+seq1_reg_rev_rtf+r')',
                                                     r'\cf2\1\cf1', i)).replace('\n', '\line '))
        file.write('\n\nsecond: ' + seq2 + '\n\n')
        file2.write('\line \line second: ' + seq2 + '\line \line ')
        for i in second:
            file.write(i)
            file2.write(re.sub(r'('+seq2_reg_rtf+r')',
                               r'\cf2\1\cf1', re.sub(r'('+seq2_reg_rev_rtf+r')',
                                                     r'\cf2\1\cf1', i)).replace('\n', '\line '))
        file.write('\n\ntogether: \n\n')
        file2.write('\line \line together: \line \line ')
        for i in together:
            file.write(i)
            tog_first = re.sub(r'('+seq1_reg_rev_rtf+r')', r'\cf2\1\cf1', i)
            tog_second = (re.sub(r'('+seq2_reg_rev_rtf+r')', r'\cf2\1\cf1', tog_first))
            tog_third = (re.sub(r'('+seq2_reg_rtf+r')', r'\cf2\1\cf1', tog_second))
            file2.write(re.sub(r'('+seq1_reg_rtf+r')', r'\cf2\1\cf1', tog_third).replace('\n', '\line '))
    file.close()
    file2.write('}')
    file2.close()

    mirnas.close()
    return [filename, filename2], status

"""parser parses pre-miRNA and mature miRNAs sequences from miRBase."""

from Bio import SeqIO
from parser_loop_localization import create_loc
import re
import RNA


def fold(sequence_to_fold):
    """
    fold folds pre-miRNA molecules and returns terminal loop.

    input: RNA sequence
    output: trimmed RNA sequence
    """
    p = re.compile('\(\.+\)')
    m = p.search(RNA.fold(sequence_to_fold)[0])
    return sequence_to_fold[m.start()-1:m.end()+1]


def parse_sh(file_name):
    """
    parse_sh parses pre-miRNA sequences for mmu ath and hsa.

    sequences from miRBase hairpin.fa file
    input: file_name with pre-miRNA sequences
    generates six files
    output: 0
    """

    fasta_sequences = SeqIO.parse(open(file_name), 'fasta')
    with open('miR_loops_mmu', 'w') as output_l_mmu, \
            open('miR_sh_mmu', 'w') as output_sh_mmu, \
            open('miR_loops_ath', 'w') as output_l_ath, \
            open('miR_sh_ath', 'w') as output_sh_ath, \
            open('miR_loops_hsa', 'w') as output_l_hsa, \
            open('miR_sh_hsa', 'w') as output_sh_hsa:
        for fasta in fasta_sequences:
            if fasta.id[:3] == 'mmu':
                output_l_mmu.write(fasta.id + ' ' +
                                   fold(str(fasta.seq)) + '\n')
                output_sh_mmu.write(fasta.id + ' ' +
                                    str(fasta.seq) + '\n')
            if fasta.id[:3] == 'ath':
                output_l_ath.write(fasta.id + ' ' +
                                   fold(str(fasta.seq)) + '\n')
                output_sh_ath.write(fasta.id + ' ' +
                                    str(fasta.seq) + '\n')
            if fasta.id[:3] == 'hsa':
                output_l_hsa.write(fasta.id + ' ' +
                                   fold(str(fasta.seq)) + '\n')
                output_sh_hsa.write(fasta.id + ' ' +
                                str(fasta.seq) + '\n')
    return 0


def parse_mature(file_name):
    """
    parse_mature parses mature miRNA sequences for mmu ath and hsa.

    sequences from miRBase mature.fa file
    input: file_name with mature miRNA sequences
    generates three files
    output: 0
    """
    fasta_sequences = SeqIO.parse(open(file_name), 'fasta')
    with open('miR_mature_mmu', 'w') as output_mature_mmu, \
            open('miR_mature_ath', 'w') as output_mature_ath, \
            open('miR_mature_hsa', 'w') as output_mature_hsa:
        for fasta in fasta_sequences:
            if fasta.id[:3] == 'mmu':
                output_mature_mmu.write(fasta.id + ' ' +
                                        str(fasta.seq) + '\n')
            if fasta.id[:3] == 'ath':
                output_mature_ath.write(fasta.id + ' ' +
                                        str(fasta.seq) + '\n')
            if fasta.id[:3] == 'hsa':
                output_mature_hsa.write(fasta.id + ' ' +
                                        str(fasta.seq) + '\n')
    return 0


if __name__ == "__main__":
    parse_sh('hairpin.fa')
    parse_mature('mature.fa')
    create_loc()

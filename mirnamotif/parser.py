"""parser parses pre-miRNA and mature miRNAs sequences from miRBase."""

from Bio import SeqIO
from paser_loop_localization import create_loc


def parse_sh(file_name):
    """
    parse_sh parses pre-miRNA sequences for mmu ath and hsa.

    sequences from miRBase hairpin.fa file
    input: file_name with pre-miRNA sequences
    generates six files
    output: 0
    """

    fasta_sequences = SeqIO.parse(open(file_name), 'fasta')
    with open('miR_sh_mmu', 'w') as output_sh_mmu, \
            open('miR_sh_ath', 'w') as output_sh_ath, \
            open('miR_sh_hsa', 'w') as output_sh_hsa:
        for fasta in fasta_sequences:
            if fasta.id[:3] == 'mmu':
                output_sh_mmu.write(fasta.id + ' ' +
                                    str(fasta.seq) + '\n')
            if fasta.id[:3] == 'ath':
                output_sh_ath.write(fasta.id + ' ' +
                                    str(fasta.seq) + '\n')
            if fasta.id[:3] == 'hsa':
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

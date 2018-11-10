

def file_to_fasta(filename):
    input_f = open(filename)
    output = open(filename[1:]+'.fa', 'w+')
    for line in input_f.readlines():
        if len(line.split()) > 1:
            output.write('>'+line.split()[0]+'\n')
            output.write(line.split()[1]+'\n')
    input_f.close()
    output.close()


def concat_files():
    output1 = open('./miR_loops.fa', 'w+')
    output2 = open('./miR_sh.fa', 'w+')
    output3 = open('./miR_mature.fa', 'w+')
    output4 = open('./miR_linking.fa', 'w+')
    l1 = open('./miR_loops_ath.fa').read()
    l2 = open('./miR_loops_hsa.fa').read()
    l3 = open('./miR_loops_mmu.fa').read()
    s1 = open('./miR_sh_ath.fa').read()
    s2 = open('./miR_sh_hsa.fa').read()
    s3 = open('./miR_sh_mmu.fa').read()
    m1 = open('./miR_mature_ath.fa').read()
    m2 = open('./miR_mature_hsa.fa').read()
    m3 = open('./miR_mature_mmu.fa').read()
    ll1 = open('./miR_linking_ath.fa').read()
    ll2 = open('./miR_linking_hsa.fa').read()
    ll3 = open('./miR_linking_mmu.fa').read()
    output1.write(l1+l2+l3)
    output2.write(s1+s2+s3)
    output3.write(m1 + m2 + m3)
    output4.write(ll1 + ll2 + ll3)
    output1.close()
    output2.close()
    output3.close()
    output4.close()


file_to_fasta('../miR_loops_ath')
file_to_fasta('../miR_loops_hsa')
file_to_fasta('../miR_loops_mmu')
file_to_fasta('../miR_sh_ath')
file_to_fasta('../miR_sh_hsa')
file_to_fasta('../miR_sh_mmu')
file_to_fasta('../miR_mature_ath')
file_to_fasta('../miR_mature_hsa')
file_to_fasta('../miR_mature_mmu')
file_to_fasta('../miR_linking_ath')
file_to_fasta('../miR_linking_hsa')
file_to_fasta('../miR_linking_mmu')
concat_files()

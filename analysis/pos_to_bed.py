import sys

def rc(seq):
    '''This function returns a reverse complement
    of a DNA or RNA strand'''

    seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()

    # reverse strand
    seq = seq[::-1]
    return seq

def get_positions(primer_set:set, rev_set:set, k:int, genome:str):
    # seq, name_list, seq_list, = '', [], []
    # for line in open(genome):
    #     line=line.strip()
    #     if '>' in line:
    #         name_list.append(line.lstrip(">").split(" ")[0])
    #         if len(seq)>0:
    #             seq_list.append(seq)
    #             seq=''
    #     else:
    #         seq=seq+line.upper()
    # seq_list.append(seq)
    # genome_list = [[name, seq_list[name_number]] for name_number, name in enumerate(name_list)]

    kmer_dict = {}
    rev_dict = {}
    for chrom, seq in genome_list:
        kmer_dict[chrom] = {}
        rev_dict[chrom] = {}
        for index in range(len(seq)-k+1):
            kmer = seq[index:index+k]
            if kmer in primer_set:
                if kmer not in kmer_dict[chrom]:
                    kmer_dict[chrom][kmer] = []
                kmer_dict[chrom][kmer].append(index)
            if kmer in rev_set:
                if kmer not in rev_dict[chrom]:
                    rev_dict[chrom][kmer] = []
                rev_dict[chrom][kmer].append(index)
    return kmer_dict, rev_dict

if __name__ == "__main__":
    genome = "/Users/kalebzuckerman/Desktop/paper_pushing/GCA_000001405.29_GRCh38.p14_genomic.fna"
    seq, name_list, seq_list, = '', [], []
    line_count = 0
    for line in open(genome, "r"):
        line_count += 1
        if line_count % 100000 == 0: print(line_count)
        line=line.strip()
        if '>' in line:
            name_list.append(line.lstrip(">").split(" ")[0])
            if len(seq)>0:
                seq_list.append(seq)
                seq=''
        else:
            seq=seq+line.upper()
    seq_list.append(seq)
    global genome_list
    genome_list = [[name, seq_list[name_number]] for name_number, name in enumerate(name_list)]
    print("opened genome")

    coat = ['ATATATATATATTA', 'ATATATATTTATTT', 'TTATATTTTTTTTT', 'CCTATATTAGTTA', 'ATATAGGTCTTA', 'AATATATTTTTTTTT', 'TCATTATTATTATC', 'TATATATTTATTTA', 'TATATTTTTTTTAT', 'ATTATTTATATTATT', 'TATATTATATTTTTT', 'ATTATTCATATTATT', 'TATATTTTTTTTTTA']
    sanger = ["ATATATATATA","TATATATATATT","TATATATATAA","TAATATATAT","TATATATATTT","ATTATTATTAT","TAATAATAATA","AAAAAAAAAAAA","AATAATAATAA","TATTATATAT"]
    soap = ["TTTTTGTGAT","TCATATAAAAT","TGAACGAT","CATAACTAA","ATTAAGGAT","AAGGTCATA","TCACTAATA","AAAGAATACA","AACATAAGT","ACATAAGTA","GGGATATAT","TTTTTGTGAT","GATTTAAAAT","TTAGAAAATT","GGTATTAAA","TCTACAATA","TATTGTAGA","AAAATCTTTA","TAAAGATTTT","ATGTAGTAA","TTTAAATCAA","TTGATTTAAA","TATTTTCATG","TAGTAACAT","TCTAATACA"]
    all_p = [coat, sanger, soap]
    # files = [("pf_coat_fwd.bed", "pf_coat_rev.bed"), ("pf_sang_fwd.bed", "pf_sang_rev.bed"), ("pf_soap_fwd.bed", "pf_soap_rev.bed")]
    files = [("hs_coat_fwd.bed", "hs_coat_rev.bed"), ("hs_sang_fwd.bed", "hs_sang_rev.bed"), ("hs_soap_fwd.bed", "hs_soap_rev.bed")]
    for i in range(3):
        primes = all_p[i]
        names = {primer: primes.index(primer) for primer in primes}
        # with open(sys.argv[1], 'w') as outf, open(sys.argv[2], 'w') as outr:
        dir = "/Users/kalebzuckerman/Desktop/paper_pushing/igv_vis/"
        with open(dir + files[i][0], 'w') as outf, open(dir + files[i][1], 'w') as outr:
            for k in range(6, 23):
                primers = set([primer for primer in primes if len(primer) == k])
                revs = set([rc(primer) for primer in primes if len(primer) == k])
                if len(primers) > 0:
                    print(f"{k}: {primers}")
                    fwd, rev = get_positions(primers, revs, k, "/Users/kalebzuckerman/Desktop/paper_pushing/3D7_unmasked_genome.fa")
                    for chr in fwd:
                        for primer in fwd[chr]:
                            for pos in fwd[chr][primer]:
                                outf.write(f"{chr}\t{pos}\t{pos + len(primer)}\tp{names[primer]}\n")
                        for reverse in rev[chr]:
                            prim = rc(reverse)
                            for pos in rev[chr][reverse]:
                                outr.write(f"{chr}\t{pos}\t{pos + len(reverse)}\tp{names[prim]}\n")


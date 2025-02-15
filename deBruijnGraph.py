from utils import *
from DNAdroplet import DNADroplet
import khmer
from chromosomeECC import ChromosomeECC
from chromosomeChunk import ChromosomeChunk
import crc16pure
import crc8

class DeBruijnGraph:
    def __init__(self, kmer_len=31):
        self.kmer_len = kmer_len
        self.counts = khmer.Counttable(kmer_len, 1e8, 4)

        self.p_kmer = {}
        self.p_cov = 0
        # self.cRule = {}
        self.next_bases = {}
        self.pre_bases = {}

        self.paths = []
        self.crc_paths = []

        self.index_byte_num = 4
        self.path_len = 0
        # self.dataEncodingDnaLength = 150
        self.ratio_tolerance = 20

        self.set_max_kmer_num = False
        self.max_path_num = 100000
        self.max_num_of_repeat_kmer = 4
        self.min_kmer_cov = 5

        self.ecc_mode = 'crc16' #crc16
        self.ecc_byte_num = 2


        # Just for test only, will be removed later
        self.primerF = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
        self.primerE = 'CTGACACTGATGCATCCG'  # complement seq of P2


    # 20200807 modified to fix a bug with quanlity scores
    def count_file(self, file):
        self.counts.consume_seqfile(file)

    def set_kmer_len(self, km_len, size=1e8, table_num=4):
        self.kmer_len = km_len
        self.counts = khmer.Counttable(km_len, size, table_num)

    def set_ecc_mode(self, mode):
        if mode == 'crc8':
            self.ecc_byte_num = 1
            self.ecc_mode = mode
        elif mode == 'crc16':
            self.ecc_byte_num = 2
            self.ecc_mode = mode

    def max_strand_num(self, block_seqs, min_km_cov=0):
        #maximal strand reconstruction number
        msn = 0
        kms = []
        c_path = True
        for seq in block_seqs:
            kms = kmers_of_str(seq, self.kmer_len)
            c_path = True
            for km in kms:
                if not self.counts.get(km) > min_km_cov:
                    c_path = False
            if c_path:
                msn = msn + 1
        return msn


    def assemble_chunk(self, chunk_info=ChromosomeChunk(0)):
        anchor_dna = chunk_info.anchor_dna()
        block_len = chunk_info.block_len
        chunk_size = chunk_info.chunk_size
        # chunk_len = block_len * chunk_size
        self.obtained_paths = []
        self.paths = [anchor_dna]
        # print(anchor_dna)
        # self.path_len = len(anchor_dna) + chunk_len


        # print(len(chunk_info.block_crc_codes))

        block_seqs = [] #
        assembled_block = 0

        pos = len(anchor_dna)
        i = 0
        assembled_block_path_num = 1
        while i < chunk_size and assembled_block_path_num > 0:
            assembled_block_path_num = self.assemble_block(chunk_info.block_crc_codes[i], pos, block_len)
            if assembled_block_path_num > 0:
                assembled_block = assembled_block + 1
                if len(self.paths) == 1:
                    while len(block_seqs) < assembled_block:
                        start = -( assembled_block - len(block_seqs)    ) * block_len
                        end   = len(self.paths[0]) -( assembled_block - len(block_seqs) -1 ) * block_len
                        # print("Start:" + str(start) + " End:" + str(end))
                        block_seqs.append(self.paths[0][start:end])
                else:
                    print(self.paths)
            pos = 0
            i = i + 1

        # cons_seq = self.consensus_seqs(self.paths)
        # block_seqs = self.split_chunk_seq_into_block_seqs(cons_seq, block_len, chunk_size)

        # print(len(block_seqs))

        rev_block_seqs = []
        # print("chunk_info.block_crc_codes: ", end='')

        #If some blocks cannot be recovered
        if len(block_seqs) < chunk_size:
            rev_assembled_block_path_num = 1
            next_anchor_dna = chunk_info.next_anchor_dna()
            if len(next_anchor_dna) > 0:
                self.paths = [next_anchor_dna]
                  #
                rev_block_seqs = []
                rev_assembled_block = 0
                i = 0
                pos=0
                while i < chunk_size - len(block_seqs) -len(rev_block_seqs) and rev_assembled_block_path_num > 0:
                    rev_assembled_block_path_num = self.assemble_block_rev(chunk_info.block_crc_codes[chunk_size-i-1], pos, block_len)
                    if rev_assembled_block_path_num > 0:
                        rev_assembled_block = rev_assembled_block + 1
                        if len(self.paths) == 1:
                            while len(rev_block_seqs) < rev_assembled_block:
                                start = (rev_assembled_block-len(rev_block_seqs)-1) * block_len #-( assembled_block - len(block_seqs)    ) * block_len
                                end   = (rev_assembled_block-len(rev_block_seqs)) * block_len #len(self.paths[0]) -( assembled_block - len(block_seqs) -1 ) * block_len
                                # print("Start:" + str(start) + " End:" + str(end))
                                rev_block_seqs.insert(0, self.paths[0][start:end])
                    i = i + 1

        if len(rev_block_seqs) > 0:
            print("********** Reverse found blocks: " + str(len(rev_block_seqs)) +  "*************************")

        #adding AAAAAAAAAAAA
        while len(block_seqs) + len(rev_block_seqs) < chunk_size:
            #block_seqs.append("A" * block_len)
            block_seqs.append("")

        return block_seqs + rev_block_seqs



    def assemble_block(self, block_crc=bytes(1), start_pos=0, block_len=500):
        # block_len = self.block_len
        # chunk_size = self.chunk_size
        # print("assemble_block() function block length: " + str(block_len))
        i = start_pos
        while i < block_len:
            tmpPaths = []
            for path in self.paths:
                # print(path)
                path_kms = kmers_sta_str(path, self.kmer_len)

                self.p_kmer = path[len(path) - self.kmer_len:len(path)]
                # print("p_kmer: " + self.p_kmer)
                self.p_cov = self.get_kmer_cov(self.p_kmer)

                maxScore = self.score_next_bases()
                for base in self.next_bases.keys():
                    if self.next_bases[base] > self.min_kmer_cov and ratio(self.next_bases[base], self.p_cov)  <= self.ratio_tolerance and len(tmpPaths) <= self.max_path_num: # and self.kmer_freq_in_path(self.p_kmer[1:] + base, path) <= self.max_num_of_repeat_kmer:
                        #20250208 Lifu Song
                        next_kmer = self.p_kmer[1:self.kmer_len] + base[0:1]
                        if next_kmer in path_kms.keys():
                            if path_kms[next_kmer] < self.max_num_of_repeat_kmer:
                                tmpPaths.append(path + base)
                            else:
                                print("Too many repeated k-mers")
                        else:
                            tmpPaths.append(path + base)

                    # print(i, end='\t')
                    # print(self.p_kmer, end='\t')
                    # print(self.p_kmer[1:] + base, end='\t')
                    # print("Path number:", end='\t')
                    # print(len(tmpPaths))
                    # print(tmpPaths)

            if len(tmpPaths) > self.max_path_num:
                # print(tmpPaths)
                tmpPaths = []
                print("Too many paths found!")

            if len(tmpPaths) > 0:
                self.paths = tmpPaths
            else:
                return False
            i = i + 1

        # check blocks one by one
        tmpPaths = []
        for path in self.paths:
            if len(path) >= block_len:
                b_seq = path[-block_len:]
                if self.ecc_mode == 'crc16':
                    if crc16pure.crc16xmodem(DNAToBytes(b_seq)) == int.from_bytes(block_crc, 'big', signed=False):
                        tmpPaths.append(path)
                else:
                    if self.ecc_mode == 'crc8':
                        crc8_obj = crc8.crc8()
                        crc8_obj.update(DNAToBytes(b_seq))
                        if crc8_obj.digest() == block_crc: #int.from_bytes(block_crc, 'big', signed=False):
                            tmpPaths.append(path)
                # else:
                #     print(path)
                #     print("CRC16 fails")
        if len(tmpPaths) > 0:
            # print("\nTotal paths found in block:" + str(len(self.paths)), end='\t')
            self.paths = tmpPaths
            # print("CRC passed paths:" + str(len(self.paths)), end='\n')

        return len(tmpPaths)


    def assemble_block_rev(self, block_crc=bytes(1), start_pos=0, block_len=500):
        # block_len = self.block_len
        # chunk_size = self.chunk_size
        # print("assemble_block_rev() function block length: " + str(block_len))
        i = start_pos
        while i < block_len:
            tmpPaths = []
            for path in self.paths:
                # print(path)
                self.p_kmer = path[:self.kmer_len]  #path[len(path) - self.kmer_len:len(path)]

                # print("Hello")
                # print(path)
                # print(self.p_kmer)
                # print(len(self.p_kmer))

                self.p_cov = self.get_kmer_cov(self.p_kmer)

                maxScore = self.score_pre_bases()
                for block in self.pre_bases.keys():
                    if self.pre_bases[block] > self.min_kmer_cov and ratio(self.pre_bases[block], self.p_cov)  <= self.ratio_tolerance and len(tmpPaths) <= self.max_path_num:
                        tmpPaths.append(block + path)

            if len(tmpPaths) > self.max_path_num:
                tmpPaths = []
                print("Too many paths found!")

            if len(tmpPaths) > 0:
                self.paths = tmpPaths
            else:
                return False
            i = i + 1

        # check blocks one by one
        tmpPaths = []
        for path in self.paths:
            if len(path) >= block_len:
                b_seq = path[:block_len]
                if crc16pure.crc16xmodem(DNAToBytes(b_seq)) == int.from_bytes(block_crc, 'big', signed=False):
                    tmpPaths.append(path)
                # else:
                #     print(path)
                #     print("CRC16 fails")

        if len(tmpPaths) > 0:
            # print("\nTotal paths found in block:" + str(len(self.paths)), end='\t')
            self.paths = tmpPaths
            # print("CRC passed paths:" + str(len(self.paths)), end='\n')

        return len(tmpPaths)


    def score_next_bases(self):
        self.next_bases = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        max_score = 0
        for m in self.next_bases.keys():
            self.next_bases[m] = self.score_next_base(m)
            if self.next_bases[m] > max_score:
                max_score = self.next_bases[m]
        return max_score


    def score_next_base(self, aBase):
        p_kmer = self.p_kmer
        next_kmer = p_kmer[1:len(p_kmer)] + aBase[0:1]
        return self.get_kmer_cov(next_kmer)


    def score_pre_bases(self):
        self.pre_bases = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        max_score = 0
        for m in self.pre_bases.keys():
            self.pre_bases[m] = self.score_pre_base(m)
            if self.pre_bases[m] > max_score:
                max_score = self.pre_bases[m]
        return max_score

    def score_pre_base(self, aBase):
        p_kmer = self.p_kmer
        pre_kmer = aBase[0:1] + p_kmer[0:-1]
        return self.get_kmer_cov(pre_kmer)

    def get_kmer_cov(self, kmer):
        if len(kmer) == self.kmer_len:
            a = self.counts.get(kmer)
            b = self.counts.get(DNA_rev_complement(kmer))
            # print("kmer cov a: ", end="")
            # print(a)
            # print(" kmer cov b: ", end="")
            # print(b)
            if kmer != DNA_rev_complement(kmer):
                return a + b
            else:
                return a
        else:
            return 0








    # def split_chunk_seq_into_block_seqs(self, chunk_seq, block_len=500, chunk_size=6):
    #     block_num = math.ceil(len(chunk_seq)/block_len)
    #     block_seqs = []
    #     i = 0
    #     for i in range(0, block_num):
    #         block_seqs.append(chunk_seq[i * block_len])
    #
    #     while len(block_seqs) < chunk_size:
    #         block_seqs.append("")
    #        # block_seqs.insert(0, "")


    # def split_chunk_seq_into_block_seqs_rev(self, chunk_seq, block_len=500, chunk_size=6):
    #     block_num = math.ceil(len(chunk_seq)/block_len)
    #     block_seqs = []
    #     i = 0
    #     for i in range(0, block_num):
    #         block_seqs.insert(0, chunk_seq[-(i+1)*block_len:-i*block_len])
    #
    #     while len(block_seqs) < chunk_size:
    #         #block_seqs.append("")
    #         block_seqs.insert(0, "")
    #
    # def consensus_seqs(self, seqs):
    #     con_seq = ""
    #     seq_len = len(seqs[0])
    #     seq_num = len(seqs)
    #     i = 0
    #     for i in range(0, seq_len):
    #         j = 0
    #         bs = seqs[0][i:i+1]
    #         for j in range(1, seq_num):
    #             if not bs == seqs[j][i:i+1]:
    #                 return con_seq
    #         con_seq = con_seq + bs
    #     return con_seq
    #
    # def consensus_seqs_rev(self, seqs):
    #     con_seq = ""
    #     seq_len = len(seqs[0])
    #     seq_num = len(seqs)
    #     i = 0
    #     for i in range(0, seq_len):
    #         j = 0
    #         bs = seqs[0][-i-1:-i]
    #         for j in range(1, seq_num):
    #             if not bs == seqs[j][-i-1:-i]:
    #                 return con_seq
    #         con_seq = bs + con_seq
    #     return con_seq
    #


    def kmer_freq_in_path(self, kmer, path):
        kmers = {}
        kmer_len = len(kmer)
        if len(path) >= kmer_len:
            i = 0
            kmstr = ''
            while i <= len(path) - kmer_len:
                kmstr = path[i:i + kmer_len]
                if kmstr in kmers.keys():
                    kmers[kmstr] = kmers[kmstr] + 1
                else:
                    if DNA_rev_complement(kmstr) in kmers.keys():
                        kmers[DNA_rev_complement(kmstr)] = kmers[DNA_rev_complement(kmstr)] + 1
                    else:
                        kmers[kmstr] = 1
                i = i + 1
        if kmer in kmers.keys():
            return kmers[kmer]
        else:
            if DNA_rev_complement(kmer) in kmers.keys():
                return kmers[DNA_rev_complement(kmer)]
            else:
                return 0
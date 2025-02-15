from utils import *
from test_utils import *
from chromosomeChunk import ChromosomeChunk
import numpy as np
import array
from reedsolo import RSCodec, ReedSolomonError
import zlib
import crc16pure
from hashlib import md5
import crc8
import re
from chromosomeECC import ChromosomeECC

class GenomeECC:
    def __init__(self, block_len=300, chunk_size=10, redundant_block_rate=0.04):
        self.block_len = block_len
        self.chunk_size = chunk_size
        self.chunk_len = self.block_len * self.chunk_size
        self.redundant_block_rate = redundant_block_rate

        # Shared variants
        self.genome_seqs = {}
        self.chr_ECC_objs = {}
        self.chr_ecc_arr = []

    def read_genome(self, file):
        self.genome_seqs = self.read_fasta(file)
        self.chr_num = len(self.genome_seqs)

    def read_chromosome(self, file):
        chro_seqs = self.read_fasta(file)
        for seq_head in chro_seqs:
            self.add_chromosome(seq_head, chro_seqs[seq_head])

    def add_chromosome(self, chro_head, chro_seq):
        self.genome_seqs[chro_head] = chro_seq

    def genChromesomeECCs(self):
        for chr_head in self.genome_seqs:
            self.chr_ECC_objs[chr_head] = self.genChromesomeECC(chr_head, self.genome_seqs[chr_head])


    def genChromesomeECC(self, head_info, seq):
        circular_chr = False
        if 'circular' in head_info:
            circular_chr = True
        aChrECC = ChromosomeECC(self.block_len, self.chunk_size, self.redundant_block_rate)
        aChrECC.circular_DNA = circular_chr
        aChrECC.set_chr_seq(seq)
        return aChrECC

    def genome_ecc(self):
        # chr_ecc_len = 0
        # chr_ecc_bytes = b''
        # self.DNA_seq_len.to_bytes(4, "big", signed=False)
        ecc_bytes_arr = [self.chr_num.to_bytes(1, "big", signed=False)]
        for head_info in self.chr_ECC_objs:
            chr_ecc_bytes = self.chr_ECC_objs[head_info].get_all_ecc_bytes()
            chr_ecc_len = len(chr_ecc_bytes)

            ecc_bytes_arr.append(chr_ecc_len.to_bytes(4, "big", signed=False))
            ecc_bytes_arr.append(chr_ecc_bytes)
        self.chr_ecc_arr = ecc_bytes_arr
        return ecc_bytes_arr

    def get_all_genome_ecc(self):
        return b''.join(self.genome_ecc())

    def write_ecc_to_file(self, filename):
        bytes_data = self.get_all_genome_ecc()
        with open(filename, "wb") as file:
            file.write(bytes_data)

    def read_ecc_file(self, filename):
        with open(filename, "rb") as file:
            file_bytes = file.read()
        return file_bytes

    def genomeEcc_to_chromosomeECC_arr(self, ecc_bytes):
        chr_num = int.from_bytes(ecc_bytes[:1], 'big', signed=False)
        print("Total number of ECC bytes: ", end='')
        print(len(ecc_bytes))
        print("Chrosome num: ", end='')
        print(chr_num)


        chr_ecc_arr = []
        curr_pos = 1
        for i in range(0, chr_num):
            b_len = int.from_bytes(ecc_bytes[curr_pos:curr_pos+4], 'big', signed=False)
            print("Current chrosome ECC length: ", end='')
            print(b_len)
            print("Current pos: ", end='')
            print(curr_pos)

            chr_ecc_arr.append(ecc_bytes[curr_pos+4:curr_pos+4+b_len])
            curr_pos = curr_pos + 4 + b_len
        self.chr_ecc_arr = chr_ecc_arr

    def read_fasta(self, file):
        f = open(file, "r")
        seq_hash = {}
        seqs = []
        matchLineA = re.compile('^(>)')

        line = f.readline()
        seq = ''
        seq_name = ''
        while line.strip():
            if not matchLineA.match(line):
                seq = seq + line.strip()
            else:
                if seq != '':
                    seqs.append(seq)
                    seq_hash[seq_name] = seq
                    seq = ''
                seq_name = self.extract_id(line)
            line = f.readline()
        f.close()
        if seq != '':
            seqs.append(seq)
            seq_hash[seq_name] = seq
        return seq_hash

    def extract_id(self, input_string):
        # Define a regular expression pattern to match the ID
        pattern = r'\|([^\|]+)\|'
        # Use re.search to find the first match in the input string
        match = re.search(pattern, input_string)
        # Check if a match is found
        if match:
            # Group 1 of the match contains the ID
            id_value = match.group(1)
            return id_value
        else:
            # Return None if no match is found
            return None
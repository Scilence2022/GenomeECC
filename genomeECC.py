import logging
import re

from chromosomeECC import ChromosomeECC

logger = logging.getLogger(__name__)


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
        # Full '>' header per chromosome id, kept so circular-chromosome detection can
        # look at the description rather than the bare accession.
        self.chr_headers = {}

    @property
    def chr_num(self):
        # Was a plain attribute set only by read_genome(), so genome_ecc() raised
        # AttributeError after read_chromosome()/add_chromosome().
        return len(self.genome_seqs)

    def read_genome(self, file):
        self.genome_seqs = {}
        self.chr_headers = {}
        self.genome_seqs = self.read_fasta(file)

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
        # head_info is the chromosome id; the 'circular' keyword lives in the full FASTA
        # description, so consult that when we have it. Guarded against None, which used
        # to raise TypeError on any header without a |...| accession.
        full_header = self.chr_headers.get(head_info) or (head_info or '')
        circular_chr = 'circular' in full_header.lower()
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
        logger.debug("Total ECC bytes: %d, chromosomes: %d", len(ecc_bytes), chr_num)

        chr_ecc_arr = []
        curr_pos = 1
        for i in range(0, chr_num):
            b_len = int.from_bytes(ecc_bytes[curr_pos:curr_pos+4], 'big', signed=False)
            logger.debug("chromosome %d: ECC length %d at pos %d", i, b_len, curr_pos)
            chr_ecc_arr.append(ecc_bytes[curr_pos+4:curr_pos+4+b_len])
            curr_pos = curr_pos + 4 + b_len
        self.chr_ecc_arr = chr_ecc_arr

    def read_fasta(self, file):
        # Rewritten: the previous "while line.strip()" loop exited at the first blank
        # line, silently truncating multi-record genomes, and built each sequence with
        # repeated string concatenation (quadratic on multi-Mbp chromosomes).
        seq_hash = {}
        seq_name = None
        header = ''
        parts = []

        def flush():
            # Merges rather than resets, so read_chromosome() can accumulate records
            # across several files.
            if seq_name is not None and parts:
                seq_hash[seq_name] = ''.join(parts)
                self.chr_headers[seq_name] = header

        with open(file, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    flush()
                    parts = []
                    header = line
                    seq_name = self.extract_id(line)
                    if seq_name is None or seq_name in seq_hash:
                        seq_name = line[1:].split()[0] if len(line) > 1 else 'seq_%d' % len(seq_hash)
                    while seq_name in seq_hash:
                        seq_name += '_dup'
                else:
                    parts.append(line)
            flush()
        return seq_hash

    def extract_id(self, input_string):
        # Prefer a |...| accession field, else fall back to the first token of the
        # header. Returning None for anything without pipes meant every standard NCBI
        # record ('>NC_000913.3 Escherichia coli ...') collapsed onto the same None key,
        # silently keeping only the last chromosome of a multi-chromosome genome.
        match = re.search(r'\|([^\|]+)\|', input_string)
        if match:
            return match.group(1)
        token = input_string.lstrip('>').split()
        return token[0] if token else None
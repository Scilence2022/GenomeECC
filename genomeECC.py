"""Compatibility API for whole-genome legacy ECC archives."""

from pathlib import Path

from chromosomeECC import ChromosomeECC
from ecc_format import unpack_legacy_genome, write_atomic
from fasta import extract_id, extract_topology, read_fasta


class GenomeECC:
    def __init__(self, block_len=300, chunk_size=10, redundant_block_rate=0.04):
        self.block_len = block_len
        self.chunk_size = chunk_size
        self.chunk_len = block_len * chunk_size
        self.redundant_block_rate = redundant_block_rate
        self.genome_seqs = {}
        self.chr_ECC_objs = {}
        self.chr_ecc_arr = []
        self.chr_headers = {}

    @property
    def chr_num(self):
        return len(self.genome_seqs)

    def read_genome(self, file):
        records = list(read_fasta(file))
        self.genome_seqs = {record.name: record.sequence for record in records}
        self.chr_headers = {record.name: f">{record.title}" for record in records}
        self.chr_ECC_objs = {}
        self.chr_ecc_arr = []

    def read_chromosome(self, file):
        records = list(read_fasta(file))
        for record in records:
            if record.name in self.genome_seqs:
                raise ValueError(f"duplicate FASTA identifier {record.name!r}")
        for record in records:
            self.add_chromosome(record.name, record.sequence)
            self.chr_headers[record.name] = f">{record.title}"

    def add_chromosome(self, chro_head, chro_seq):
        self.genome_seqs[chro_head] = chro_seq
        self.chr_headers.pop(chro_head, None)
        self.chr_ECC_objs.pop(chro_head, None)
        self.chr_ecc_arr = []

    def genChromesomeECCs(self):
        self.chr_ECC_objs = {
            head: self.genChromesomeECC(head, sequence) for head, sequence in self.genome_seqs.items()
        }
        self.chr_ecc_arr = []

    def genChromesomeECC(self, head_info, seq):
        title = self.chr_headers.get(head_info) or str(head_info or '')
        chromosome = ChromosomeECC(self.block_len, self.chunk_size, self.redundant_block_rate)
        chromosome.circular_DNA = extract_topology(title) == 'circular'
        chromosome.set_chr_seq(seq)
        return chromosome

    def genome_ecc(self):
        if not 1 <= self.chr_num <= 255:
            raise ValueError("legacy genome ECC supports 1..255 records")
        if self.chr_ECC_objs.keys() != self.genome_seqs.keys() or any(
            self.chr_ECC_objs[head].DNA_seq != seq for head, seq in self.genome_seqs.items()
        ):
            raise ValueError("generate ECC for the current genome before writing")
        parts = [bytes([self.chr_num])]
        for head in self.genome_seqs:
            chromosome = self.chr_ECC_objs[head].get_all_ecc_bytes()
            parts.extend([len(chromosome).to_bytes(4, "big"), chromosome])
        self.chr_ecc_arr = parts
        return parts

    def get_all_genome_ecc(self):
        return b''.join(self.genome_ecc())

    def write_ecc_to_file(self, filename):
        write_atomic(filename, self.get_all_genome_ecc())

    def read_ecc_file(self, filename):
        return Path(filename).read_bytes()

    def genomeEcc_to_chromosomeECC_arr(self, ecc_bytes):
        self.chr_ecc_arr = [payload for payload, _ in unpack_legacy_genome(ecc_bytes)]

    def read_fasta(self, file):
        records = list(read_fasta(file))
        self.chr_headers.update({record.name: f">{record.title}" for record in records})
        return {record.name: record.sequence for record in records}

    def extract_id(self, input_string):
        return extract_id(input_string)

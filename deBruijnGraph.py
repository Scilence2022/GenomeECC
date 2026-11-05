"""Bounded, bidirectional block assembly from canonical k-mer counts."""

import logging

from chromosomeChunk import ChromosomeChunk, block_crc
from chromosomeECC import ChromosomeECC
from ecc_format import CRC_WIDTHS
from utils import DNAToBytes, rev_seq


logger = logging.getLogger(__name__)


def _khmer_table(kmer_len, size, tables):
    try:
        import khmer
    except ImportError as exc:
        raise ImportError(
            "read assembly requires khmer 3.0.0a3; install it from Bioconda "
            "(see README), or recover from --blocks"
        ) from exc
    return khmer.Counttable(kmer_len, size, tables)


class DeBruijnGraph:
    def __init__(self, kmer_len=31, *, counter=None):
        if not 1 <= kmer_len <= 32:
            raise ValueError("k-mer length must be in 1..32")
        self.kmer_len = kmer_len
        self.counts = _khmer_table(kmer_len, 1e8, 4) if counter is None else counter
        self.paths = []
        self.ratio_tolerance = 20
        self.max_path_num = 1000
        self.max_num_of_repeat_kmer = 4
        self.min_kmer_cov = 5
        self.set_ecc_mode("crc16")

    def count_file(self, file):
        self.counts.consume_seqfile(file)

    def set_kmer_len(self, km_len, size=1e8, table_num=4):
        if not 1 <= km_len <= 32:
            raise ValueError("k-mer length must be in 1..32")
        self.counts = _khmer_table(km_len, size, table_num)
        self.kmer_len = km_len

    def set_ecc_mode(self, mode):
        if mode not in CRC_WIDTHS:
            raise ValueError(f"unknown block_ecc_mode {mode!r}")
        self.ecc_mode = mode
        self.ecc_byte_num = CRC_WIDTHS[mode]

    def get_kmer_cov(self, kmer):
        if len(kmer) != self.kmer_len:
            return 0
        return self.counts.get(min(kmer, rev_seq(kmer))) or 0

    def max_strand_num(self, block_seqs, min_km_cov=0):
        k = self.kmer_len
        return sum(
            len(seq) >= k and all(
                self.get_kmer_cov(seq[i:i + k]) > min_km_cov for i in range(len(seq) - k + 1)
            )
            for seq in block_seqs
        )

    def _extend(self, steps, reverse):
        k = self.kmer_len
        for _ in range(steps):
            candidates = []
            for path in self.paths:
                if len(path) < k:
                    continue
                previous = path[:k] if reverse else path[-k:]
                coverage = self.get_kmer_cov(previous)
                for base in "ACGT":
                    kmer = base + previous[:-1] if reverse else previous[1:] + base
                    next_coverage = self.get_kmer_cov(kmer)
                    if next_coverage < self.min_kmer_cov:
                        continue
                    if coverage and max(coverage, next_coverage) > self.ratio_tolerance * min(
                        coverage, next_coverage
                    ):
                        continue
                    if self.kmer_freq_in_path(kmer, path) >= self.max_num_of_repeat_kmer:
                        continue
                    if len(candidates) == self.max_path_num:
                        logger.debug("assembly exceeded the path budget of %d", self.max_path_num)
                        self.paths = []
                        return
                    candidates.append(base + path if reverse else path + base)
            self.paths = candidates
            if not candidates:
                return

    def _assemble_block(self, checksum, steps, length, padded_len, reverse, mode):
        self._extend(steps, reverse)
        self.paths = [
            path for path in self.paths if len(path) >= length and block_crc(
                DNAToBytes((path[:length] if reverse else path[-length:]).ljust(padded_len, 'A')),
                mode,
            ) == checksum
        ]
        return len(self.paths)

    def assemble_block(self, block_crc=bytes(1), start_pos=0, block_len=500):
        return self._assemble_block(block_crc, block_len - start_pos, block_len, block_len,
                                    False, self.ecc_mode)

    def assemble_block_rev(self, block_crc=bytes(1), start_pos=0, block_len=500):
        return self._assemble_block(block_crc, block_len - start_pos, block_len, block_len,
                                    True, self.ecc_mode)

    def assemble_chunk(self, chunk_info=None):
        if chunk_info is None:
            chunk_info = ChromosomeChunk(0)
        count, block_len = chunk_info.chunk_size, chunk_info.block_len
        length = chunk_info.valid_bases
        if length is None:
            length = count * block_len
        if not (count - 1) * block_len < length <= count * block_len:
            raise ValueError("invalid chunk geometry")
        if len(chunk_info.block_crc_codes) != count:
            raise ValueError("chunk checksum count does not match its geometry")
        mode = chunk_info.ecc_mode if chunk_info.valid_bases is not None else self.ecc_mode
        lengths = [min(block_len, length - i * block_len) for i in range(count)]
        anchor = chunk_info.anchor_dna()[:lengths[0]]
        self.paths = [anchor]
        blocks = []
        for i, width in enumerate(lengths):
            steps = width - len(anchor) if i == 0 else width
            if not self._assemble_block(chunk_info.block_crc_codes[i], steps, width, block_len,
                                        False, mode):
                break
            if len(self.paths) == 1:
                blocks = [
                    self.paths[0][j * block_len:j * block_len + lengths[j]] for j in range(i + 1)
                ]

        reverse_blocks = []
        next_anchor = chunk_info.next_anchor_dna()
        self.paths = [next_anchor] if len(next_anchor) >= self.kmer_len else []
        for i in range(count - 1, len(blocks) - 1, -1):
            width = lengths[i]
            if not self._assemble_block(chunk_info.block_crc_codes[i], width, width, block_len,
                                        True, mode):
                break
            if len(self.paths) == 1:
                reverse_blocks = [
                    self.paths[0][(j - i) * block_len:(j - i) * block_len + lengths[j]]
                    for j in range(i, count)
                ]
        missing = count - len(blocks) - len(reverse_blocks)
        if missing:
            logger.debug("chunk %d has %d unassembled blocks", chunk_info.chunk_index, missing)
        return blocks + [''] * missing + reverse_blocks

    def kmer_freq_in_path(self, kmer, path):
        occurrences = 0
        position = path.find(kmer)
        while position != -1 and occurrences < self.max_num_of_repeat_kmer:
            occurrences += 1
            position = path.find(kmer, position + 1)
        return occurrences


def assemble_candidates(archive, graph):
    """Return global block positions and unambiguous DNA candidates for each record."""
    result = {}
    for record in archive.records:
        chunks = []
        for group in record.groups:
            chromosome = ChromosomeECC()
            chromosome.circular_DNA = False
            chromosome.reset_with_unpacked_ecc(group.header)
            chunks.extend(chromosome.rebuild_chunk_objs_with_ecc_codes(group.header["all_chunk_ecc_bytes"]))
        for i, chunk in enumerate(chunks):
            if i + 1 < len(chunks) or record.topology == "circular":
                following = chunks[(i + 1) % len(chunks)]
                chunk.next_chunk_anchor = following.chunk_anchor
                chunk.next_anchor_valid_bases = following.valid_bases
        result[record.name] = candidates = {}
        offset = 0
        for chunk in chunks:
            for i, sequence in enumerate(graph.assemble_chunk(chunk), offset):
                if sequence:
                    candidates[i] = sequence
            offset += chunk.chunk_size
    return result

import crc8
import crc16pure

from ecc_format import CRC_WIDTHS
from utils import DNAToBytes, bytesToDNA


def block_crc(payload, mode):
    if mode == "crc8":
        return crc8.crc8(payload).digest()
    if mode == "crc16":
        return crc16pure.crc16xmodem(payload).to_bytes(2, "big")
    raise ValueError(f"unknown block_ecc_mode {mode!r}")


class ChromosomeChunk:
    def __init__(self, chunk_index=1, chunk_anchor=b'', crc_codes=b'', block_len=500, chunk_size=6):
        self.chunk_index = chunk_index
        self.chunk_anchor = chunk_anchor
        self.next_chunk_anchor = b''
        self.next_anchor_valid_bases = None
        self.chunk_size = chunk_size
        self.set_block_len(block_len)
        self.valid_bases = None
        self.block_crc_codes = []
        self.crc_codes = crc_codes
        self.set_ecc_mode('crc16')
        self.circular = True

    def set_ecc_mode(self, mode):
        if mode not in CRC_WIDTHS:
            raise ValueError(f"unknown block_ecc_mode {mode!r}")
        self.ecc_byte_num = CRC_WIDTHS[mode]
        self.ecc_mode = mode

    def anchor_dna(self):
        dna = bytesToDNA(self.chunk_anchor)
        return dna if self.valid_bases is None else dna[:self.valid_bases]

    def next_anchor_dna(self):
        return bytesToDNA(self.next_chunk_anchor)[:self.next_anchor_valid_bases]

    def set_block_len(self, block_len):
        self.block_len = block_len
        self.block_size = block_len // 4
        self.chunk_len = block_len * self.chunk_size

    def check_block_crc(self, block_seq):
        data = DNAToBytes(block_seq)
        if len(data) < self.ecc_byte_num:
            return False
        return block_crc(data[:-self.ecc_byte_num], self.ecc_mode) == data[-self.ecc_byte_num:]

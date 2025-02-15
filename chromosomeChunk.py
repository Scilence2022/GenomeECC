from utils import *
from DNAdroplet import DNADroplet
import khmer
import numpy as np
import array
from reedsolo import RSCodec, ReedSolomonError
# rsc = RSCodec(10)  # 10 ecc symbols
import zlib
import crc16pure
import crc8
from hashlib import md5



class ChromosomeChunk:
    def __init__(self, chunk_index=1, chunk_anchor=b'', crc_codes=b'', block_len=500, chunk_size=6):

        self.chunk_index = chunk_index
        self.chunk_anchor = chunk_anchor
        self.next_chunk_anchor = b''
        self.block_len = block_len
        self.block_size = int(self.block_len/4)
        self.chunk_size = chunk_size
        self.chunk_len = self.block_len * self.chunk_size

        self.block_crc_codes = []
        self.crc_codes = crc_codes
        self.ecc_mode = 'crc16' #crc16
        self.ecc_byte_num = 2
        self.circular = True

    def set_ecc_mode(self, mode):
        if mode == 'crc8':
            self.ecc_byte_num = 1
            self.ecc_mode = mode
        elif mode == 'crc16':
            self.ecc_byte_num = 2
            self.ecc_mode = mode


    def anchor_dna(self):
        return bytesToDNA(self.chunk_anchor)

    def next_anchor_dna(self):
        return bytesToDNA(self.next_chunk_anchor)

    def set_block_len(self, block_len):
        self.block_len = block_len
        self.block_size = int(self.block_len / 4)
        self.chunk_len = self.block_len * self.chunk_size

    def check_block_crc(self, block_seq):
        block_bytes = bytesToDNA(block_seq)
        if self.ecc_mode == 'crc16':
            if crc16pure.crc16xmodem(block_bytes[:-self.ecc_byte_num]) == block_bytes[-self.ecc_byte_num:]:
                return True
            else:
                return False
        else:
            if self.ecc_mode == 'crc8':
                crc8_obj = crc8.crc8()
                crc8_obj.update(block_bytes[:-self.ecc_byte_num])

                if crc8_obj.digest() == block_bytes[-self.ecc_byte_num:]:
                    return True
                else:
                    return False
        # def get_block_ecc(self, block)

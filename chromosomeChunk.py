from utils import *
import crc16pure
import crc8


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
        # block_seq is DNA whose trailing ecc_byte_num bytes carry the block's CRC.
        # This used to call bytesToDNA (the wrong direction, raising TypeError) and
        # then compare an int against bytes, which is never equal.
        block_bytes = DNAToBytes(block_seq)
        payload = block_bytes[:-self.ecc_byte_num]
        stored = block_bytes[-self.ecc_byte_num:]
        if self.ecc_mode == 'crc16':
            expected = crc16pure.crc16xmodem(payload).to_bytes(self.ecc_byte_num, 'big', signed=False)
        elif self.ecc_mode == 'crc8':
            crc8_obj = crc8.crc8()
            crc8_obj.update(payload)
            expected = crc8_obj.digest()
        else:
            raise ValueError("unknown ecc_mode %r" % (self.ecc_mode,))
        return expected == stored
        # def get_block_ecc(self, block)

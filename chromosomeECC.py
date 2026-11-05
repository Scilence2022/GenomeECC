"""Reed–Solomon protection for one independently encoded sequence or group."""

from hashlib import md5
import math
from pathlib import Path
from threading import RLock
import zlib

import numpy as np
from reedsolo import RSCodec, ReedSolomonError

from chromosomeChunk import ChromosomeChunk, block_crc
from ecc_format import (
    CRC_MODES, CRC_WIDTHS, MAX_CODEWORD, RS_POLYNOMIAL, ceil_div,
    unpack_chromosome, validate_codeword, validate_geometry, validate_rate, write_atomic,
)
from utils import DNAToBytes, bytesToDNA


# reedsolo 1.x uses module-global finite-field tables during construction and coding.
_RS_LOCK = RLock()


class ChromosomeECC:
    nsize = MAX_CODEWORD
    colum_byte_num_per_block = 2

    def __init__(self, block_len=500, chunk_size=6, redundant_block_rate=0.03,
                 block_ecc_mode='crc8', chunk_anchor_byte_num=8):
        self.block_len = block_len
        self.block_size = block_len // 4
        self.chunk_size = chunk_size
        self.chunk_len = block_len * chunk_size
        self.redundant_block_rate = redundant_block_rate
        self.chunk_anchor_byte_num = chunk_anchor_byte_num
        self.circular_DNA = True

        self.DNA_seq_len = 0
        self.block_num = 0
        self.chunk_num = 0
        self.redundant_block_num = 0
        self.ecc = None
        self._clear_sequence()
        self.set_block_ecc_mode(block_ecc_mode)

    def _clear_sequence(self):
        self.DNA_seq = ''
        self.DNA_seq_bytes = b''
        self.DNA_seq_bytes_md5 = None
        self.block_crc_codes = []
        self.chunk_crc_codes = []
        self._loaded_anchors = []
        self.block_seqs = []
        self.block_bytes = []
        self.colum_bytes = []
        self.colum_rs_codes = []

    def _geometry(self, length, block_len, chunk_size):
        validate_geometry(block_len, chunk_size, self.chunk_anchor_byte_num)
        validate_rate(self.redundant_block_rate)
        if type(length) is not int or not 0 < length < 2**32:
            raise ValueError("chromosome length must be in 1..4294967295 bases")
        blocks = ceil_div(length, block_len)
        parity = math.ceil(self.redundant_block_rate * blocks)
        validate_codeword(blocks, parity)
        return blocks, parity

    def _configure(self, length, blocks, parity):
        self.DNA_seq_len = length
        self.block_num = blocks
        self.chunk_num = ceil_div(blocks, self.chunk_size)
        self._set_parity(parity)

    def _set_parity(self, parity):
        self.redundant_block_num = parity
        self.colum_rs_codes = []
        self.ecc = None

    def _codec(self):
        if self.ecc is None:
            self.ecc = RSCodec(
                self.redundant_block_num, nsize=self.nsize, prim=RS_POLYNOMIAL,
                generator=2, fcr=0, c_exp=16,
            )
        return self.ecc

    def set_chr_seq(self, seq):
        blocks, parity = self._geometry(len(seq), self.block_len, self.chunk_size)
        if self.block_ecc_mode not in CRC_WIDTHS:
            raise ValueError(f"unknown block_ecc_mode {self.block_ecc_mode!r}")
        packed = DNAToBytes(self.round_DNA(seq))
        self._clear_sequence()
        self._configure(len(seq), blocks, parity)
        self.DNA_seq = seq
        self.DNA_seq_bytes = packed
        self.chunk_crc_codes = [None] * self.chunk_num
        self.get_block_crc(0)

    def set_chr_size(self, g_size):
        blocks, parity = self._geometry(g_size, self.block_len, self.chunk_size)
        self._clear_sequence()
        self._configure(g_size, blocks, parity)
        self.chunk_crc_codes = [None] * self.chunk_num

    def _set_geometry(self, block_len, chunk_size):
        validate_geometry(block_len, chunk_size, self.chunk_anchor_byte_num)
        if self.DNA_seq_bytes_md5 is not None:
            raise ValueError("cannot change the geometry of a loaded recovery context")
        if self.DNA_seq_len:
            self._geometry(self.DNA_seq_len, block_len, chunk_size)
        same_blocks = self.block_len == block_len
        self.block_len = block_len
        self.block_size = block_len // 4
        self.chunk_size = chunk_size
        self.chunk_len = block_len * chunk_size
        if self.DNA_seq:
            if same_blocks:
                self.chunk_num = ceil_div(self.block_num, chunk_size)
                self.chunk_crc_codes = [None] * self.chunk_num
            else:
                self.set_chr_seq(self.DNA_seq)
        elif self.DNA_seq_len:
            self.set_chr_size(self.DNA_seq_len)

    def set_block_len(self, block_len):
        self._set_geometry(block_len, self.chunk_size)

    def set_chunk_size(self, chunk_size):
        self._set_geometry(self.block_len, chunk_size)

    def set_block_ecc_mode(self, mode):
        if mode not in CRC_WIDTHS:
            raise ValueError(f"unknown block_ecc_mode {mode!r}")
        if self.DNA_seq_bytes_md5 is not None and mode != self.block_ecc_mode:
            raise ValueError("cannot change the CRC mode of a loaded recovery context")
        self.block_ecc_mode = mode
        self.block_ecc_byte_num = CRC_WIDTHS[mode]
        if self.DNA_seq:
            self.block_crc_codes = []

    def set_redundant_block_rate(self, redundant_block_rate):
        validate_rate(redundant_block_rate)
        self.set_redundant_block_num(math.ceil(redundant_block_rate * self.block_num))
        self.redundant_block_rate = redundant_block_rate

    def set_redundant_block_num(self, redundant_block_num):
        validate_codeword(self.block_num, redundant_block_num)
        if self.DNA_seq_bytes_md5 is not None and redundant_block_num != self.redundant_block_num:
            raise ValueError("cannot change the parity of a loaded recovery context")
        self._set_parity(redundant_block_num)
        # This setter has historically made the codec immediately available to callers.
        with _RS_LOCK:
            self._codec()

    @property
    def chunk_anchores(self):
        if self.DNA_seq:
            return [self.get_chunk_anchor(i) for i in range(self.chunk_num)]
        return self._loaded_anchors

    def get_chunk_anchor(self, chunk_id):
        if not self.DNA_seq:
            return self._loaded_anchors[chunk_id]
        return self.chunk(chunk_id)[:self.chunk_anchor_byte_num].ljust(self.chunk_anchor_byte_num, b'\0')

    def get_chunk_index(self, chunk_id):
        return self.get_chunk_anchor(chunk_id) or None

    def get_chunk_block_crc_codes(self, chunk_id):
        start = chunk_id * self.chunk_size
        end = min(start + self.chunk_size, self.block_num)
        codes = [self.get_block_crc(i) for i in range(start, end)]
        return codes + [bytes(self.block_ecc_byte_num)] * (self.chunk_size - len(codes))

    def round_DNA(self, dnaseq='', group_size=4):
        return dnaseq + 'A' * (-len(dnaseq) % group_size)

    def chunk(self, chunk_id):
        start = self.chunk_size * chunk_id * self.block_size
        return self.DNA_seq_bytes[start:start + self.chunk_size * self.block_size]

    def chunk_seq(self, chunk_id):
        return bytesToDNA(self.chunk(chunk_id))

    def block_seqs_to_bytes(self):
        if len(self.block_seqs) != self.block_num:
            raise ValueError(f"expected {self.block_num} blocks, got {len(self.block_seqs)}")
        blocks = []
        for i, sequence in enumerate(self.block_seqs):
            actual_len = min(self.block_len, self.DNA_seq_len - i * self.block_len)
            if len(sequence) not in (actual_len, self.block_len):
                raise ValueError(f"block {i}: expected {actual_len} bases")
            blocks.append(DNAToBytes(sequence[:actual_len].ljust(self.block_len, 'A')))
        self.block_bytes = blocks
        self.colum_bytes = []
        self.colum_rs_codes = []

    def block(self, block_id):
        if not 0 <= block_id < self.block_num:
            return b''
        if not self.block_bytes:
            if not self.DNA_seq:
                raise ValueError("no sequence blocks loaded")
            packed = self.DNA_seq_bytes
            size = self.block_size
            self.block_bytes = [
                packed[start:start + size].ljust(size, b'\0')
                for start in range(0, len(packed), size)
            ]
        return self.block_bytes[block_id]

    def colum(self, colum_id):
        if not self.colum_bytes:
            self.block_bytes_to_colum_bytes()
        return self.colum_bytes[colum_id]

    def colum_num(self):
        return ceil_div(self.block_size, self.colum_byte_num_per_block)

    def _padded_block_bytes(self, block_id):
        return self.block_bytes[block_id].ljust(self.colum_num() * self.colum_byte_num_per_block, b'\0')

    def block_bytes_to_colum_bytes(self):
        self.block(0)
        if len(self.block_bytes) != self.block_num or any(
            len(block) != self.block_size for block in self.block_bytes
        ):
            raise ValueError("blocks do not match the configured geometry")
        matrix = np.frombuffer(b''.join(self.block_bytes), dtype=np.uint8)
        matrix = matrix.reshape(self.block_num, self.block_size)
        if self.block_size % 2:
            matrix = np.pad(matrix, ((0, 0), (0, 1)))
        self.colum_bytes = [matrix[:, start:start + 2].tobytes() for start in range(0, matrix.shape[1], 2)]
        self.colum_rs_codes = []

    def restore_colum_bytes_rs(self, rs_bytes, erasure_positions=()):
        validate_codeword(self.block_num, self.redundant_block_num)
        positions = list(erasure_positions)
        if any(type(pos) is not int or not 0 <= pos < self.block_num for pos in positions):
            raise ValueError("erasure positions must be valid block indices")
        if len(set(positions)) != len(positions):
            raise ValueError("duplicate erasure position")
        if len(positions) > self.redundant_block_num:
            raise ReedSolomonError(
                f"too many erasures: {len(positions)} for {self.redundant_block_num} parity symbols"
            )
        parity = self._split_parity(rs_bytes)
        if not self.colum_bytes:
            self.block_bytes_to_colum_bytes()
        if len(self.colum_bytes) != self.colum_num() or any(
            len(col) != 2 * self.block_num for col in self.colum_bytes
        ):
            raise ValueError("columns do not match the configured geometry")
        if not self.redundant_block_num:
            return 0
        restored = []
        corrected_symbols = 0
        with _RS_LOCK:
            codec = self._codec()
            for column, parity_column in zip(self.colum_bytes, parity):
                symbols = np.frombuffer(column + parity_column, dtype='>u2').tolist()
                decoded, _, errata = codec.decode(symbols, erase_pos=positions)
                restored.append(np.asarray(decoded, dtype='>u2').tobytes())
                corrected_symbols += len(errata)
        self.colum_bytes = restored
        self.colum_rs_codes = parity
        return corrected_symbols

    def _split_parity(self, rs_bytes):
        width = self.colum_byte_num_per_block * self.redundant_block_num
        if len(rs_bytes) != self.colum_num() * width:
            raise ValueError("Reed–Solomon parity length does not match the configured geometry")
        return [rs_bytes[i * width:(i + 1) * width] for i in range(self.colum_num())]

    def reset_rs_bytes_with_rs_codes(self, rs_bytes):
        self.colum_rs_codes = self._split_parity(rs_bytes)

    def colum_bytes_to_block_bytes(self):
        columns = np.frombuffer(b''.join(self.colum_bytes), dtype=np.uint8)
        matrix = columns.reshape(self.colum_num(), self.block_num, 2).transpose(1, 0, 2)
        matrix = matrix.reshape(self.block_num, -1)
        self.block_bytes = [row[:self.block_size].tobytes() for row in matrix]

    def block_seq(self, block_id):
        return bytesToDNA(self.block(block_id))

    def get_block_seqs(self):
        return [self.block_seq(i) for i in range(self.block_num)]

    def write_block_seqs_to_fasta(self, file):
        fasta = ''.join(f">block-{i}\n{self.block_seq(i)}\n" for i in range(self.block_num))
        write_atomic(file, fasta.encode('ascii'))

    def chunk_block_seq(self, chunk_id, block_index):
        return self.block_seq(chunk_id * self.chunk_size + block_index)

    def get_chr_md5(self):
        if self.DNA_seq_bytes_md5 is not None:
            return self.DNA_seq_bytes_md5
        return md5(self.DNA_seq_bytes).digest()

    def get_chunk_crc(self, chunk_id):
        if self.chunk_crc_codes[chunk_id] is None:
            self.chunk_crc_codes[chunk_id] = zlib.crc32(self.chunk(chunk_id)).to_bytes(4, 'big')
        return self.chunk_crc_codes[chunk_id]

    def get_block_crc(self, block_id):
        if not self.block_crc_codes:
            self.block_crc_codes = [
                block_crc(self.block(i), self.block_ecc_mode) for i in range(self.block_num)
            ]
        return self.block_crc_codes[block_id]

    def get_colum_rs(self, colum_id):
        if not self.colum_rs_codes:
            self.__colum_rs()
        return self.colum_rs_codes[colum_id]

    def __colum_rs(self):
        validate_codeword(self.block_num, self.redundant_block_num)
        if not self.redundant_block_num:
            self.colum_rs_codes = [b''] * self.colum_num()
            return
        if not self.colum_bytes:
            self.block_bytes_to_colum_bytes()
        parity = []
        with _RS_LOCK:
            codec = self._codec()
            for column in self.colum_bytes:
                symbols = np.frombuffer(column, dtype='>u2').tolist()
                encoded = codec.encode(symbols)
                parity.append(np.asarray(encoded[self.block_num:], dtype='>u2').tobytes())
        self.colum_rs_codes = parity

    def get_chunk_ecc(self, chunk_id):
        return self.get_chunk_anchor(chunk_id) + b''.join(self.get_chunk_block_crc_codes(chunk_id))

    def pack_anchor_crc_byte_num(self):
        return bytes([(self.chunk_anchor_byte_num << 4) | self.block_ecc_byte_num])

    def pack_block_len_chunk_size(self):
        if not 0 <= self.block_len < 1024:
            raise ValueError(f"block_len must fit in 10 bits, got {self.block_len!r}")
        if not 0 <= self.chunk_size < 64:
            raise ValueError(f"chunk_size must fit in 6 bits, got {self.chunk_size!r}")
        return ((self.block_len << 6) | self.chunk_size).to_bytes(2, 'big')

    def unpack_block_len_chunk_size(self, packed_bytes=b''):
        packed = int.from_bytes(packed_bytes, 'big')
        return [packed >> 6, packed & 63]

    def chr_ecc(self):
        if not self.DNA_seq:
            raise ValueError("load a reference sequence before encoding")
        validate_geometry(self.block_len, self.chunk_size, self.chunk_anchor_byte_num)
        if (
            self.block_size != self.block_len // 4
            or self.block_num != ceil_div(self.DNA_seq_len, self.block_len)
            or self.chunk_len != self.block_len * self.chunk_size
            or self.chunk_num != ceil_div(self.block_num, self.chunk_size)
            or self.block_ecc_byte_num != CRC_WIDTHS.get(self.block_ecc_mode)
        ):
            raise ValueError("chromosome configuration changed without its setter")
        validate_codeword(self.block_num, self.redundant_block_num)
        chunks = [self.get_chunk_ecc(i) for i in range(self.chunk_num)]
        parity = [self.get_colum_rs(i) for i in range(self.colum_num())]
        chunk_width = self.chunk_anchor_byte_num + self.chunk_size * self.block_ecc_byte_num
        if any(len(chunk) != chunk_width for chunk in chunks) or any(
            len(column) != 2 * self.redundant_block_num for column in parity
        ):
            raise ValueError("stored checksums or parity do not match the chromosome geometry")
        return [
            self.DNA_seq_len.to_bytes(4, 'big'),
            self.get_chr_md5(),
            self.pack_block_len_chunk_size(),
            self.pack_anchor_crc_byte_num(),
            self.redundant_block_num.to_bytes(2, 'big'),
            *chunks,
            *parity,
        ]

    def unpack_ecc_bytes(self, ecc_bytes):
        return unpack_chromosome(ecc_bytes)

    def reset_with_unpacked_ecc(self, unpacked_ecc):
        self._clear_sequence()
        self.DNA_seq_len = 0
        self.block_len = unpacked_ecc['block_len']
        self.chunk_size = unpacked_ecc['chunk_size']
        self.chunk_len = self.block_len * self.chunk_size
        self.block_size = unpacked_ecc['block_size']
        self.chunk_anchor_byte_num = unpacked_ecc['chunk_anchor_byte_num']
        self.set_block_ecc_mode(CRC_MODES[unpacked_ecc['block_crc_byte_num']])
        self._configure(
            unpacked_ecc['genome_length'], unpacked_ecc['block_num'],
            unpacked_ecc['redundant_block_num'],
        )
        self.DNA_seq_bytes_md5 = unpacked_ecc['genome_md5']
        self.chunk_crc_codes = [None] * self.chunk_num
        if "circular" in unpacked_ecc:
            self.circular_DNA = unpacked_ecc["circular"]

    def rebuild_chunk_objs_with_ecc_codes(self, chunk_ecc_bytes):
        anchor_width, crc_width = self.chunk_anchor_byte_num, self.block_ecc_byte_num
        width = anchor_width + crc_width * self.chunk_size
        if len(chunk_ecc_bytes) != width * self.chunk_num:
            raise ValueError("chunk metadata length does not match the configured geometry")
        chunks, checksums = [], []
        for i in range(self.chunk_num):
            start = i * width
            anchor = chunk_ecc_bytes[start:start + anchor_width]
            remaining = self.DNA_seq_len - i * self.chunk_len
            valid_bases = min(self.chunk_len, remaining)
            valid_blocks = ceil_div(valid_bases, self.block_len)
            chunk = ChromosomeChunk(i, anchor, block_len=self.block_len, chunk_size=valid_blocks)
            chunk.valid_bases = valid_bases
            chunk.set_ecc_mode(self.block_ecc_mode)
            chunk.circular = self.circular_DNA
            start += anchor_width
            chunk.block_crc_codes = [
                chunk_ecc_bytes[start + j * crc_width:start + (j + 1) * crc_width]
                for j in range(valid_blocks)
            ]
            if chunks:
                chunks[-1].next_chunk_anchor = anchor
                chunks[-1].next_anchor_valid_bases = valid_bases
            chunks.append(chunk)
            checksums.extend(chunk.block_crc_codes)
        if self.circular_DNA and chunks:
            chunks[-1].next_chunk_anchor = chunks[0].chunk_anchor
            chunks[-1].next_anchor_valid_bases = chunks[0].valid_bases
        self.block_crc_codes = checksums
        self._loaded_anchors = [chunk.chunk_anchor for chunk in chunks]
        return chunks

    def get_all_ecc_bytes(self):
        return b''.join(self.chr_ecc())

    def write_ecc_to_file(self, filename):
        write_atomic(filename, self.get_all_ecc_bytes())

    def read_ecc_file(self, filename):
        return Path(filename).read_bytes()

    def assembled_chr_seq(self):
        self.block(0)
        return bytesToDNA(b''.join(self.block_bytes))[:self.DNA_seq_len]

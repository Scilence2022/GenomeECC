from utils import *
from test_utils import *
# from DNAdroplet import DNADroplet
# import khmer
from chromosomeChunk import ChromosomeChunk
import numpy as np
import array
from reedsolo import RSCodec, ReedSolomonError
import zlib
import crc16pure
from hashlib import md5
import crc8

class ChromosomeECC:
    def __init__(self, block_len=500, chunk_size=6, redundant_block_rate=0.03, block_ecc_mode='crc8', chunk_anchor_byte_num=8):
        self.block_len = block_len
        self.block_size = int(self.block_len/4)
        self.chunk_size = chunk_size
        self.chunk_len = self.block_len * self.chunk_size
        self.redundant_block_rate = redundant_block_rate
        self.chunk_anchor_byte_num = chunk_anchor_byte_num

        self.DNA_seq = ''
        self.DNA_seq_len = 0
        self.circular_DNA = True
        self.block_num = 0
        self.chunk_num = 0
        self.colum_byte_num_per_block = 2
        self.redundant_block_num = 0

        self.DNA_seq_bytes = b''
        # block_ecc_byte_num must follow block_ecc_mode: crc8 digests are 1 byte. It was
        # hardcoded to 2 while the default mode is 'crc8', so the packed header declared
        # a 2-byte CRC width for 1-byte codes and every unpack misaligned.
        self.block_ecc_mode = block_ecc_mode
        self.block_ecc_byte_num = 2
        self.set_block_ecc_mode(block_ecc_mode)
        self.block_crc_codes = [] #np.array(bytes(2))
        # object dtype: np.array(bytes(4)) built an |S4 array that silently truncated
        # any CRC32 ending in a NUL byte.
        self.chunk_crc_codes = np.empty(0, dtype=object)
        self.chunk_anchores = []


        self.DNA_seq_bytes_md5 = ""
        #md5(some_bytes).hexdigest()
        self.ecc = RSCodec()
        self.nsize = int(math.pow(2,16)) - 1 # 65535

        self.block_bytes = []
        self.colum_bytes = []
        self.colum_rs_codes = []
        self.block_seqs = []

    def block_seqs_to_bytes(self):
        self.block_bytes = []
        for seq in self.block_seqs:
            self.block_bytes.append(DNAToBytes(seq))

    def set_chr_seq(self, seq):
        self.DNA_seq = seq
        self.DNA_seq_len = len(seq)
        self.DNA_seq_bytes = DNAToBytes(self.round_DNA(self.DNA_seq))
        self.block_num = math.ceil(self.DNA_seq_len/self.block_len)
        self.chunk_num = math.ceil(self.DNA_seq_len/self.chunk_len)
        # Reset derived state first: __block_crc and __encode_chunk_anchores append, so
        # calling set_chr_seq twice used to double these lists.
        self.block_crc_codes = []
        self.chunk_anchores = []
        self.block_bytes = []
        self.colum_bytes = []
        self.colum_rs_codes = []
        self.chunk_crc_codes = np.empty(self.chunk_num, dtype=object)
        self.__block_crc()
        self.__chunk_crc()
        self.__encode_chunk_anchores()

        self.redundant_block_num = math.ceil(self.redundant_block_rate * self.block_num)
        self.ecc = RSCodec(self.redundant_block_num, nsize=self.nsize)


    def set_chr_size(self, g_size):
        # self.DNA_seq = seq
        self.DNA_seq_len = g_size
        # self.DNA_seq_bytes = DNAToBytes(self.round_DNA(self.DNA_seq))
        self.block_num = math.ceil(self.DNA_seq_len/self.block_len)
        self.chunk_num = math.ceil(self.DNA_seq_len/self.chunk_len)
        # block_crc_codes is a plain list; .resize() is a numpy method and raised
        # AttributeError here.
        self.block_crc_codes = [b''] * self.block_num
        self.chunk_crc_codes = np.empty(self.chunk_num, dtype=object)
        # self.chunk_anchores.resize(self.chunk_num)

        self.redundant_block_num = math.ceil(self.redundant_block_rate * self.block_num)
        self.ecc = RSCodec(self.redundant_block_num, nsize=self.nsize)


    def set_block_len(self, block_len):
        self.block_len = block_len
        self.block_size = int(self.block_len / 4)
        # self.chunk_size = chunk_size
        self.chunk_len = self.block_len * self.chunk_size

    def set_block_ecc_mode(self, mode):
        if mode == 'crc8':
            self.block_ecc_byte_num = 1
            self.block_ecc_mode = mode
        elif mode == 'crc16':
            self.block_ecc_byte_num = 2
            self.block_ecc_mode = mode

    # nsize must be carried through: rebuilding RSCodec without it silently dropped the
    # code from GF(2^16) back to the 255-symbol default, which cannot address the
    # 2-byte-per-block symbols this class encodes.
    def set_redundant_block_rate(self, redundant_block_rate):
        self.redundant_block_rate = redundant_block_rate
        self.redundant_block_num = math.floor(redundant_block_rate*self.block_num)
        self.ecc = RSCodec(self.redundant_block_num, nsize=self.nsize)

    def set_redundant_block_num(self, redundant_block_num):
        self.redundant_block_num = math.floor(redundant_block_num)
        self.ecc = RSCodec(self.redundant_block_num, nsize=self.nsize)


    def set_chunk_size(self, chunk_size):
        self.chunk_size = chunk_size
        ##

    def get_chunk_anchor(self, chunk_id):
        return self.chunk(chunk_id)[:self.chunk_anchor_byte_num]


    def get_chunk_index(self, chunk_id):
        if self.chunk_anchores[chunk_id]:
            return self.chunk_anchores[chunk_id]

    def get_chunk_block_crc_codes(self, chunk_id):
        start = chunk_id * self.chunk_size
        end = (chunk_id + 1) * self.chunk_size
        if end > self.block_num:
            # Pad the short final chunk with zero CRCs. np.resize tiled the existing
            # codes cyclically instead of padding, and coerced them to a fixed-width
            # bytes dtype that stripped any CRC ending in a NUL byte.
            block_crc_codes = list(self.block_crc_codes[start:self.block_num])
            filler = bytes(self.block_ecc_byte_num)
            block_crc_codes += [filler] * (self.chunk_size - len(block_crc_codes))
        else:
            block_crc_codes = list(self.block_crc_codes[start:end])

        return block_crc_codes

    def __encode_chunk_anchores(self):
        # Was hardcoded to 6 bytes while get_chunk_anchor (and the packed header) use
        # chunk_anchor_byte_num, which defaults to 8 -- the two disagreed.
        self.chunk_anchores = []
        for i in range(0, self.chunk_num):
            self.chunk_anchores.append(self.chunk(i)[0:self.chunk_anchor_byte_num])

    def round_DNA(self, dnaseq='', group_size=4):
        # % group_size on the padding width: when the length was already aligned this
        # appended a whole spurious group (4 bases), inflating DNA_seq_bytes by a byte
        # and folding that byte into the stored MD5.
        return dnaseq + ('A' * ((group_size - len(dnaseq) % group_size) % group_size))


    def chunk(self, chunk_id):
        start = self.chunk_size * chunk_id * self.block_size
        end = min(self.chunk_size * self.block_size * (chunk_id + 1), len(self.DNA_seq_bytes))
        return self.DNA_seq_bytes[start:end]

    def chunk_seq(self, chunk_id):
        return bytesToDNA(self.chunk(chunk_id))

    def block(self, block_id):
        if len(self.block_bytes) < self.block_num:
            self._gen_block_bytes()
        if block_id < self.block_num:
            return self.block_bytes[block_id]
        else:
            return b''

    def _block_byte(self, block_id):
        start = self.block_size * block_id
        end = min(self.block_size * (block_id + 1), len(self.DNA_seq_bytes))
        if end-start < self.block_size:
            return self.DNA_seq_bytes[start:end] + bytes(self.block_size - (end-start) )
        return self.DNA_seq_bytes[start:end]

    def _gen_block_bytes(self):
        self.block_bytes = []
        for i in range(0, self.block_num):
            self.block_bytes.append(self._block_byte(i))
        # while len(self.block_bytes[-1]


    def colum(self, colum_id):
        #generating block_bytes
        if len(self.colum_bytes) < self.colum_num():
            self.block_bytes_to_colum_bytes()
        return self.colum_bytes[colum_id]


    def colum_num(self):
        # ceil, not int(): with the default block_len=500 the block is 125 bytes, and
        # truncating 125/2 to 62 columns covered only 124 bytes. The last byte of every
        # block was in no column at all -- neither RS-protected nor restored on
        # reassembly -- silently losing 4 bases per block. Blocks are zero-padded up to
        # colum_num()*colum_byte_num_per_block so every column still holds exactly one
        # symbol per block, which is what keeps RS erasure positions aligned to blocks.
        return math.ceil(self.block_size / self.colum_byte_num_per_block)


    def _padded_block_bytes(self, block_id):
        b = self.block_bytes[block_id]
        pad = self.colum_num() * self.colum_byte_num_per_block - len(b)
        return b + bytes(pad) if pad > 0 else b


    def block_bytes_to_colum_bytes(self):
        if len(self.block_bytes) < self.block_num:
            self.block_bytes = []
            self._gen_block_bytes()
        n = self.colum_byte_num_per_block
        self.colum_bytes = []
        for i in range(0, self.colum_num()):
            self.colum_bytes.append(
                b''.join(self._padded_block_bytes(j)[i*n:(i+1)*n] for j in range(0, self.block_num)))

    def restore_colum_bytes_rs(self, rs_bytes, erasure_positions=[]):
        self.reset_rs_bytes_with_rs_codes(rs_bytes)
        for i in range(0, self.colum_num()):
            int_arr = bytes_to_int_array(self.colum_bytes[i] + self.colum_rs_codes[i], self.colum_byte_num_per_block)
            dec_obj = self.ecc.decode(int_arr, self.ecc.nsym, erasure_positions)
            self.colum_bytes[i] = int_array_to_bytes(dec_obj[0], self.colum_byte_num_per_block)


    def reset_rs_bytes_with_rs_codes(self, rs_bytes):
        self.colum_rs_codes = []
        rs_code_len = self.colum_byte_num_per_block * self.redundant_block_num
        for i in range(0, self.colum_num()):
            self.colum_rs_codes.append(rs_bytes[i*rs_code_len:(i+1)*rs_code_len])

    def colum_bytes_to_block_bytes(self):
        n = self.colum_byte_num_per_block
        for i in range(0, self.block_num):
            full = b''.join(self.colum_bytes[j][i*n:(i+1)*n] for j in range(0, self.colum_num()))
            self.block_bytes[i] = full[:self.block_size]   # strip the column padding

    def block_seq(self, block_id):
        return bytesToDNA(self.block(block_id))

    def get_block_seqs(self):
        sqs = []
        for i in range(0, self.block_num):
            sqs.append(self.block_seq(i))
        return sqs

    def write_block_seqs_to_fasta(self, file):
        # 2020/06/30
        out = open(file, 'tw')
        for i in range(0, self.block_num):
            out.write(">block-" + str(i) + "\n")
            out.write(self.block_seq(i))
            out.write("\n")
        out.close()

    def chunk_block_seq(self, chunk_id, block_index):
        block_id = chunk_id * self.chunk_size + block_index
        return self.block_seq(block_id)


    def __chunk_crc(self):
        for i in range(0, self.chunk_num):
            self.chunk_crc_codes[i] = zlib.crc32(self.chunk(i)).to_bytes(4, 'big', signed=False)

    def __block_crc(self):
            # Rebuilds from scratch; this used to append unconditionally, so any second
            # call (e.g. a repeated set_chr_seq) doubled the list.
            self.block_crc_codes = []
            if self.block_ecc_mode == 'crc8':
                for i in range(0, self.block_num):
                    crc8_obj = crc8.crc8()
                    crc8_obj.update(self.block(i))
                    self.block_crc_codes.append(crc8_obj.digest())
            elif self.block_ecc_mode == 'crc16':
                for i in range(0, self.block_num):
                    self.block_crc_codes.append(crc16pure.crc16xmodem(self.block(i)).to_bytes(self.block_ecc_byte_num, 'big', signed=False))
            else:
                raise ValueError("unknown block_ecc_mode %r" % (self.block_ecc_mode,))



            # self.block_crc_codes.append(crc16pure.crc16xmodem(self.block(i)).to_bytes(self.block_crc_byte_num, 'big', signed=False))
            #self.block_crc_codes[i] = crc16pure.crc16xmodem(self.block(i)).to_bytes(2, 'big', signed=False)
            # print(int.from_bytes(self.block_crc_codes[i],'big', signed=False))

    def get_chr_md5(self):
        return md5(self.DNA_seq_bytes).digest()

    # These accessors test "is it already computed?", so they check for None rather than
    # truthiness: a legitimately all-zero CRC (b'\x00') is falsy and used to force a full
    # recompute of every code on each call.
    def get_chunk_crc(self, chunk_id):
        if self.chunk_crc_codes[chunk_id] is None:
            self.__chunk_crc()
        return self.chunk_crc_codes[chunk_id]

    #def chunk_index(self, ):

    def get_block_crc(self, block_id):
        if len(self.block_crc_codes) < self.block_num:
            self.__block_crc()
        return self.block_crc_codes[block_id]


    def get_colum_rs(self, colum_id):
        if len(self.colum_rs_codes) <= colum_id:
            self.__colum_rs()
        return self.colum_rs_codes[colum_id]


    def __colum_rs(self):

        colum_int_array = []
        colum_int_array_encode = []
        colum_bytes = b''
        self.colum_rs_codes = []
        #int(self.block_size/self.colum_byte_num_per_block)
        for i in range(0, self.colum_num()):
            # print("Block CRC: ", end="")
            # print(crc16pure.crc16xmodem(self.block(i)))
            colum_int_array = bytes_to_int_array(self.colum(i), self.colum_byte_num_per_block)
            colum_int_array_encode = self.ecc.encode(colum_int_array)
            self.colum_rs_codes.append(int_array_to_bytes(colum_int_array_encode[len(colum_int_array):],  self.colum_byte_num_per_block))
            #




    def get_chunk_ecc(self, chunk_id):
        #anchor 6 bytes
        #block CRC 2 bytes per block
        all_ecc_info = b''
        all_ecc_info = all_ecc_info + self.get_chunk_anchor(chunk_id)
        for crc in self.get_chunk_block_crc_codes(chunk_id):
            all_ecc_info = all_ecc_info + crc

        return all_ecc_info

    def pack_anchor_crc_byte_num(self):
        a = self.chunk_anchor_byte_num
        b = self.block_ecc_byte_num
        ff = (a<<4) + b
        return ff.to_bytes(1, 'big', signed=False)


    def pack_block_len_chunk_size(self):
        # 16-bit field: 10 bits block_len, 6 bits chunk_size. Out-of-range values used to
        # either raise a bare OverflowError (block_len >= 1024) or silently carry into
        # the block_len bits (chunk_size >= 64, e.g. 500/64 unpacked back as 501/0).
        block_len = self.block_len
        chunk_size = self.chunk_size
        if not 0 <= block_len < 1024:
            raise ValueError(
                "block_len must fit in 10 bits (0..1023), got %r" % (block_len,))
        if not 0 <= chunk_size < 64:
            raise ValueError(
                "chunk_size must fit in 6 bits (0..63), got %r" % (chunk_size,))
        ff = (block_len<<6) + chunk_size
        return ff.to_bytes(2, 'big', signed=False)

    def unpack_block_len_chunk_size(self, packed_bytes=b''):

        ff = int.from_bytes(packed_bytes, "big", signed=False)
        block_len = ff>>6
        chunk_size = ff - (block_len<<6)
        return [block_len, chunk_size]

    def chr_ecc(self):
        # returns all the required information for recovery of a genome, especially the ECC codes.
        # The
        all_ecc_info = []

        all_ecc_info.append(self.DNA_seq_len.to_bytes(4, "big", signed=False))  # genome_size
        all_ecc_info.append(self.get_chr_md5()) # MD5 of genome bytes self.get_genome_md5()
        all_ecc_info.append(self.pack_block_len_chunk_size()) #Block_length 10 bits Chunk size 6 bits  pack_anchor_crc_byte_num
        all_ecc_info.append(self.pack_anchor_crc_byte_num())  # 4 bits each  pack_anchor_crc_byte_num
        all_ecc_info.append(self.redundant_block_num.to_bytes(2, "big", signed=False)) #RS Redundancy num
        #
        for i in range(0, self.chunk_num):
            all_ecc_info.append(self.get_chunk_ecc(i))

        for i in range(0, self.colum_num()):
            all_ecc_info.append(self.get_colum_rs(i))

        return all_ecc_info


    def unpack_ecc_bytes(self, ecc_bytes):
        #
        unpacked_ecc = {}
        unpacked_ecc['genome_length'] = int.from_bytes(ecc_bytes[:4], 'big', signed=False)  #self.DNA_seq_len
        unpacked_ecc['genome_md5'] = ecc_bytes[4:4+16]
        unpacked_ecc['block_len'], unpacked_ecc['chunk_size'] = self.unpack_block_len_chunk_size(ecc_bytes[20:22])

        unpacked_ecc['chunk_anchor_byte_num'] = int.from_bytes(ecc_bytes[22:23], 'big', signed=False)
        unpacked_ecc['chunk_anchor_byte_num'] = unpacked_ecc['chunk_anchor_byte_num']>>4
        unpacked_ecc['block_crc_byte_num'] = int.from_bytes(ecc_bytes[22:23], 'big', signed=False) - (unpacked_ecc['chunk_anchor_byte_num'] <<4)


        #caculating
        unpacked_ecc['block_size']  = int(unpacked_ecc['block_len']/4)
        unpacked_ecc['redundant_block_num'] = int.from_bytes(ecc_bytes[23:25], 'big', signed=False)
        unpacked_ecc['block_num'] = math.ceil(unpacked_ecc['genome_length']/unpacked_ecc['block_len'])
        unpacked_ecc['chunk_num'] = math.ceil(unpacked_ecc['genome_length']/(unpacked_ecc['block_len'] * unpacked_ecc['chunk_size']))

        chunk_ecc_byte_len = unpacked_ecc['chunk_anchor_byte_num'] + unpacked_ecc['block_crc_byte_num'] * unpacked_ecc['chunk_size']
        # print(chunk_ecc_byte_len)
        all_chunk_ecc_bytes = ecc_bytes[25:25+chunk_ecc_byte_len*unpacked_ecc['chunk_num']]
        all_rs_code_bytes = ecc_bytes[25+chunk_ecc_byte_len*unpacked_ecc['chunk_num']:]

        unpacked_ecc['all_chunk_ecc_bytes'] = all_chunk_ecc_bytes
        unpacked_ecc['all_rs_code_bytes'] = all_rs_code_bytes

        unpacked_ecc['chunk_objs'] = []

        return unpacked_ecc
        # chunk_ecc_byte_num =
        # ck = GenomeChunk(i)
        # ck.set_block_len(gc.block_len)
        # ck.chunk_size = gc.chunk_size
        # # ck.block_len = block_len
        # ck.block_crc_codes = gc.get_chunk_block_crc_codes(i)
        # ck.chunk_anchor = gc.chunk(i)[:6]

        # self.block_crc_codes.resize(self.block_num)
        # self.chunk_crc_codes.resize(self.chunk_num)
        # self.chunk_anchores.resize(self.chunk_num)
        # self.__block_crc()
        # self.__chunk_crc()
        # self.__encode_chunk_anchores()
        #
        # self.redundant_block_num = math.ceil(self.redundant_block_rate * self.block_num)
        # self.ecc = RSCodec(self.redundant_block_num, nsize=self.nsize)

    def reset_with_unpacked_ecc(self, unpacked_ecc):
        #
        self.DNA_seq_len = unpacked_ecc['genome_length']
        genome_md5 = unpacked_ecc['genome_md5']
        self.block_len = unpacked_ecc['block_len']
        self.chunk_size = unpacked_ecc['chunk_size']
        self.chunk_len = self.block_len * self.chunk_size

        self.block_size  = unpacked_ecc['block_size']
       # self.chunk_len =
        self.redundant_block_num = unpacked_ecc['redundant_block_num']

        self.ecc = RSCodec(self.redundant_block_num, nsize=self.nsize)
        # resetting object
        self.block_num = math.ceil(self.DNA_seq_len/self.block_len)
        self.chunk_num = math.ceil(self.DNA_seq_len/self.chunk_len)

        self.block_crc_codes = []
        # self.chunk_crc_codes.resize(self.chunk_num)
        # self.chunk_anchores.resize(self.chunk_num)

    def rebuild_chunk_objs_with_ecc_codes(self, chunk_ecc_bytes):
        #test =1
        ck_objs = []
        chunk_ecc_byte_len = self.chunk_anchor_byte_num + self.block_ecc_byte_num * self.chunk_size

        #Rebuilding genome chunks
        for i in range(0, self.chunk_num):
            one_chunk_ecc = chunk_ecc_bytes[i*chunk_ecc_byte_len:(i+1)*chunk_ecc_byte_len]
            ck = ChromosomeChunk(i)
            ck.chunk_anchor = one_chunk_ecc[:self.chunk_anchor_byte_num]
            if i > 0:
                ck_objs[i-1].next_chunk_anchor = one_chunk_ecc[:self.chunk_anchor_byte_num]

            ck.set_block_len(self.block_len)
            ck.chunk_size = self.chunk_size
            ck.block_len = self.block_len
            ck.chunk_len = self.block_len * self.chunk_size
            #setting block crc codes
            for j in range(0, self.chunk_size):
                self.block_crc_codes.append(one_chunk_ecc[self.chunk_anchor_byte_num+j*self.block_ecc_byte_num:self.chunk_anchor_byte_num + (j + 1) * self.block_ecc_byte_num])
                ck.block_crc_codes.append(one_chunk_ecc[
                                            self.chunk_anchor_byte_num + j * self.block_ecc_byte_num:self.chunk_anchor_byte_num + (
                                                        j + 1) * self.block_ecc_byte_num])
            ck_objs.append(ck)

        if self.circular_DNA:  #
            ck_objs[-1].next_chunk_anchor = ck_objs[0].chunk_anchor

        return ck_objs

    def get_all_ecc_bytes(self):
        return b''.join(self.chr_ecc())

    def write_ecc_to_file(self, filename):
        bytes_data = b''.join(self.chr_ecc())
        with open(filename, "wb") as file:
            file.write(bytes_data)


    def read_ecc_file(self, filename):
        with open(filename, "rb") as file:
            file_bytes = file.read()
        return file_bytes


    def assembled_chr_seq(self):
        all_bytes = b''

        for i in range(0, self.block_num):
            all_bytes = all_bytes + self.block_bytes[i]

        return bytesToDNA(all_bytes)[:self.DNA_seq_len]










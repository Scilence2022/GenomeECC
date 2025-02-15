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
        self.block_ecc_mode = block_ecc_mode
        self.block_ecc_byte_num = 2
        self.block_crc_codes = [] #np.array(bytes(2))
        self.chunk_crc_codes = np.array(bytes(4))  # seems to be not required at all, will be deleted later
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
        # self.block_crc_codes.resize(self.block_num)
        self.chunk_crc_codes.resize(self.chunk_num)
        # self.chunk_anchores.resize(self.chunk_num)
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
        self.block_crc_codes.resize(self.block_num)
        self.chunk_crc_codes.resize(self.chunk_num)
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

    def set_redundant_block_rate(self, redundant_block_rate):
        self.redundant_block_num = math.floor(redundant_block_rate*self.block_num)
        self.ecc = RSCodec(self.redundant_block_num)

    def set_redundant_block_num(self, redundant_block_num):
        self.redundant_block_num = math.floor(redundant_block_num)
        self.ecc = RSCodec(self.redundant_block_num)


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
        # block_crc_codes = np.array(0)
        # block_crc_codes.resize(self.chunk_size)
        if end > self.block_num:
            end = self.block_num
            block_crc_codes = self.block_crc_codes[start:]
            block_crc_codes = np.resize(block_crc_codes, (self.chunk_size))
            #block_crc_codes.resize(self.chunk_size)
            # while len(block_crc_codes) < self.block_size:
            #     block_crc_codes.append([0])
        else:
            block_crc_codes = self.block_crc_codes[start:end]
            # np.resize(block_crc_codes, (self.chunk_size))

        # print("crc start: " + str(start) + " end: " + str(end))
        # print(len(block_crc_codes))
        return block_crc_codes

    def __encode_chunk_anchores(self):
        for i in range(0, self.chunk_num):
            self.chunk_anchores.append(self.chunk(i)[0:6])

    def round_DNA(self, dnaseq='', group_size=4):
        return dnaseq + ('A' * (group_size - len(dnaseq) % group_size))


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
        if len(self.colum_bytes) < int(self.block_size/self.colum_byte_num_per_block):
            self.block_bytes_to_colum_bytes()
        return self.colum_bytes[colum_id]


    def colum_num(self):
        return int(self.block_size/self.colum_byte_num_per_block)


    def block_bytes_to_colum_bytes(self):
        if len(self.block_bytes) < self.block_num:
            self.block_bytes = []
            self._gen_block_bytes()
        self.colum_bytes = []
        for i in range(0, self.colum_num()):
            c_bytes = b''
            for j in range(0, self.block_num):
                c_bytes = c_bytes + self.block_bytes[j][i*self.colum_byte_num_per_block:(i+1)*self.colum_byte_num_per_block]
            self.colum_bytes.append(c_bytes)

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
        for i in range(0, self.block_num):
            self.block_bytes[i] = b''
            for j in range(0, self.colum_num()):
                self.block_bytes[i] = self.block_bytes[i] + self.colum_bytes[j][i*self.colum_byte_num_per_block:(i+1)*self.colum_byte_num_per_block]

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
            # print("Block CRC: ", end="")
            # print(crc16pure.crc16xmodem(self.block(i)))
            if self.block_ecc_mode == 'crc8':
                for i in range(0, self.block_num):
                    crc8_obj = crc8.crc8()
                    crc8_obj.update(self.block(i))
                    self.block_crc_codes.append(crc8_obj.digest())
            else:
                if self.block_ecc_mode == 'crc16':
                    for i in range(0, self.block_num):
                        self.block_crc_codes.append(crc16pure.crc16xmodem(self.block(i)).to_bytes(self.block_ecc_byte_num, 'big', signed=False))



            # self.block_crc_codes.append(crc16pure.crc16xmodem(self.block(i)).to_bytes(self.block_crc_byte_num, 'big', signed=False))
            #self.block_crc_codes[i] = crc16pure.crc16xmodem(self.block(i)).to_bytes(2, 'big', signed=False)
            # print(int.from_bytes(self.block_crc_codes[i],'big', signed=False))

    def get_chr_md5(self):
        return md5(self.DNA_seq_bytes).digest()

    def get_chunk_crc(self, chunk_id):
        if self.chunk_crc_codes[chunk_id]:
            return self.chunk_crc_codes[chunk_id]
        else:
            self.__chunk_crc()
            return self.chunk_crc_codes[chunk_id]

    #def chunk_index(self, ):

    def get_block_crc(self, block_id):
        if self.block_crc_codes[block_id]:
            return self.block_crc_codes[block_id]
        else:
            self.__block_crc()
            return self.block_crc_codes[block_id]


    def get_colum_rs(self, colum_id):

        if len(self.colum_rs_codes)>0 and self.colum_rs_codes[colum_id]:
            return self.colum_rs_codes[colum_id]
        else:
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
        block_len = self.block_len
        chunk_size = self.chunk_size
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










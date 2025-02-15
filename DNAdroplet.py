import crc16pure
from utils import *
from droplet import Droplet
import zlib
import json
import random

class DNADroplet(Droplet):
    def __init__(self, data=b'', index=-1, num_of_chunks=100, head_index_len=16, seed=1, crc_len=16):
        super(DNADroplet, self).__init__(data, index, num_of_chunks)
        self.seed = seed
        self.index_len = head_index_len
        self.index = index

        self.crc_len = crc_len
        self.checksum = "CRC"
        self.crc = 0
        self.allbytes = b''

        self.des = False
        self.des_data = b''
        self.sec_key = ""
        self.disorder_data_bytes = b''
        
        if data:
            self.data_len = len(data) * 4
        else:
            self.data_len = len(self.data) * 4

        self.des_data_len = 32*4
        self.data_byte_order = []

    #2020.02.19 Deal with the index-->chunks problem
    def set_index(self, index_num):
        self.index = index_num
        self.symbol_index = index_num
        self.update()

    # 2020.02.19 Deal with the index-->chunks problem
    def update(self):
        self.get_chunk_nums()

    def get_head_index_dna(self):
        aDNA = self.to_DNA()
        return aDNA[0:self.index_len]

    def get_head_index(self):
        return self.index

    def _crc(self):
        #databytes = b''
        if self.des:
            # print("DNAdroplet sec_key: " + self.sec_key)
            # Only when no encrypted data
            if not self.des_data:
                self.des_data = des_en(self.data, self.sec_key)
            data_bytes = self.index.to_bytes(int(self.index_len / 4), 'big',
                                             signed=False) + self.des_data
        else:
            data_bytes = self.index.to_bytes(int(self.index_len / 4), 'big',
                                             signed=False) + self.data

        # data_bytes = self.index.to_bytes(int(self.index_len / 4), 'big', signed=False) + self.data
        #self.crc = crc16pure.crc16xmodem(data_bytes)
        self.crc = zlib.crc32(data_bytes)


    def encry_data(self):
        self.des_data = des_en(self.data, self.sec_key)
        self.des_data_len = len(self.des_data) * 4


    def decry_data(self):
        self.data = des_de(self.des_data, self.sec_key)

    def get_crc(self):
        # if self.crc <= 0:
        self._crc()
        return self.crc.to_bytes(int(self.crc_len / 4), 'big', signed=False)


    def get_byte_order(self):
        byte_len = int(self.des_data_len/4) + int(self.crc_len/4)
        self.data_byte_order = rand_order(byte_len, self.index + int.from_bytes(self.sec_key, 'big'))
        return self.data_byte_order

    def to_DNA(self):
        self._crc()
        if self.des:
            data_bytes = self.des_data \
                       + self.crc.to_bytes(int(self.crc_len / 4), 'big', signed=False)
            self.disorder_data_bytes = sort_bytes_by_order(data_bytes, self.get_byte_order())

            allbytes = self.index.to_bytes(int(self.index_len / 4), 'big', signed=False) \
                       + self.disorder_data_bytes
        else:
            allbytes = self.index.to_bytes(int(self.index_len / 4), 'big', signed=False)\
                   + self.data \
                   + self.crc.to_bytes(int(self.crc_len/4), 'big', signed=False)
        self.allbytes = allbytes

        return bytesToDNA(allbytes)

    def set_droplet_from_DNA(self, DNA_string):
        if self.des:
            self.disorder_data_bytes = DNAToBytes(DNA_string[self.index_len:self.index_len + self.des_data_len + self.crc_len])
            data_bytes = reverse_sort_bytes_by_order(self.disorder_data_bytes, self.get_byte_order())
            self.des_data = data_bytes[:-int(self.crc_len/4)]
            # print(self.index)
            # print(self.disorder_data_bytes)
            # print(data_bytes)
            # print(self.get_crc())

            if not self.get_crc() == data_bytes[-int(self.crc_len/4):]:
                return False
            self.decry_data()
        else:
            data_bytes = DNAToBytes(DNA_string[self.index_len:self.index_len + self.data_len])
            crc_bytes = DNAToBytes(DNA_string[self.index_len + self.data_len:])
            self.data = data_bytes
            if not self.get_crc() == crc_bytes:
                return False

        self.update()
        return True



import zlib

from droplet import Droplet
from utils import DNAToBytes, bytesToDNA, des_de, des_en, rand_order, reverse_sort_bytes_by_order, sort_bytes_by_order


class DNADroplet(Droplet):
    def __init__(self, data=b'', index=-1, num_of_chunks=100, head_index_len=16, seed=1, crc_len=16):
        if type(head_index_len) is not int or head_index_len < 4 or head_index_len % 4:
            raise ValueError("droplet index width must be a positive multiple of 4 bases")
        if type(crc_len) is not int or crc_len < 16 or crc_len % 4:
            raise ValueError("CRC32 width must be a multiple of 4 and at least 16 bases")
        super().__init__(data, index, num_of_chunks)
        self.seed = seed
        self.index_len = head_index_len
        self.index = index
        self.crc_len = crc_len
        self.checksum = "CRC"
        self.crc = 0
        self.allbytes = b''
        self.des = False
        self.des_data = b''
        self.sec_key = b''
        self.disorder_data_bytes = b''
        self.data_len = len(data) * 4
        self.des_data_len = (len(data) // 8 + 1) * 32
        self.data_byte_order = []

    def set_index(self, index_num):
        self.index = index_num
        self.symbol_index = index_num
        self.update()

    def update(self):
        self.get_chunk_nums()

    def get_head_index_dna(self):
        return bytesToDNA(self._index_bytes())

    def get_head_index(self):
        return self.index

    def _index_bytes(self):
        width = self.index_len // 4
        if type(self.index) is not int or not 0 < self.index < 1 << (8 * width):
            raise ValueError("droplet index does not fit the wire header")
        return self.index.to_bytes(width, 'big')

    def _crc(self):
        if self.des and not self.des_data:
            self.encry_data()
        payload = self.des_data if self.des else self.data
        self.crc = zlib.crc32(self._index_bytes() + payload)

    def encry_data(self):
        if not isinstance(self.sec_key, bytes) or len(self.sec_key) != 8:
            raise ValueError("legacy DES requires an explicit 8-byte key")
        self.des_data = des_en(self.data, self.sec_key)
        self.des_data_len = len(self.des_data) * 4

    def decry_data(self):
        self.data = des_de(self.des_data, self.sec_key)

    def get_crc(self):
        self._crc()
        return self.crc.to_bytes(self.crc_len // 4, 'big')

    def get_byte_order(self):
        self.data_byte_order = rand_order(
            (self.des_data_len + self.crc_len) // 4, self.index + int.from_bytes(self.sec_key, 'big'),
        )
        return self.data_byte_order

    def to_DNA(self):
        index = self._index_bytes()
        crc = self.get_crc()
        if self.des:
            payload = sort_bytes_by_order(self.des_data + crc, self.get_byte_order())
            self.disorder_data_bytes = payload
        else:
            payload = self.data + crc
        self.allbytes = index + payload
        return bytesToDNA(self.allbytes)

    def set_droplet_from_DNA(self, DNA_string):
        payload_length = self.des_data_len if self.des else self.data_len
        if len(DNA_string) != self.index_len + payload_length + self.crc_len:
            return False
        try:
            encoded = DNAToBytes(DNA_string)
        except ValueError:
            return False
        index_width = self.index_len // 4
        index_bytes = encoded[:index_width]
        index = int.from_bytes(index_bytes, 'big')
        if index == 0 or (self.index > 0 and self.index != index):
            return False

        data_bytes = encoded[index_width:]
        if self.des:
            if not isinstance(self.sec_key, bytes) or len(self.sec_key) != 8:
                raise ValueError("legacy DES requires an explicit 8-byte key")
            order = rand_order(len(data_bytes), index + int.from_bytes(self.sec_key, 'big'))
            data_bytes = reverse_sort_bytes_by_order(data_bytes, order)
        crc_width = self.crc_len // 4
        payload, check = data_bytes[:-crc_width], data_bytes[-crc_width:]
        crc = zlib.crc32(index_bytes + payload)
        if crc.to_bytes(crc_width, 'big') != check:
            return False

        data = des_de(payload, self.sec_key) if self.des else payload
        if len(data) * 4 != self.data_len:
            return False
        self.data = data
        self.crc = crc
        self.allbytes = encoded
        if self.des:
            self.des_data = payload
            self.disorder_data_bytes = encoded[index_width:]
            self.data_byte_order = order
        self.set_index(index)
        return True

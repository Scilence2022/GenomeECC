from DNAdroplet import DNADroplet
from fountain import Fountain
from utils import xor


class DNAFountain(Fountain):
    def __init__(self, data=b'', chunk_size=35, index=1, seed=1,
                 droplet_head_length=16, droplet_crc_length=16):
        self.droplet_head_length = droplet_head_length
        self.droplet_crc_length = droplet_crc_length
        self.des = False
        self.sec_key = b''
        super().__init__(data, chunk_size, index, seed)

    def DNAdroplet(self):
        if not self.num_of_chunks:
            raise ValueError("empty payload has no droplets")
        self.updateIndex()
        droplet = DNADroplet(
            b'', self.index, self.num_of_chunks, self.droplet_head_length,
            self.seed, self.droplet_crc_length,
        )
        droplet.degree = self.degrees[self.index % len(self.degrees)]
        data = bytes(self.chunk_size)
        for index in droplet.get_chunk_nums():
            data = xor(data, self.chunk(index))
        droplet.data = data
        droplet.data_len = self.chunk_size * 4
        if self.des:
            droplet.des = True
            droplet.sec_key = self.sec_key
            droplet.encry_data()
        return droplet

    def fix_bytes(self):
        """Pad the backing payload for legacy callers that compare source chunks."""
        self.data += bytes(-len(self.data) % self.chunk_size)

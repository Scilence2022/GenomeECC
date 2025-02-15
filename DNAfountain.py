from DNAdroplet import DNADroplet
from fountain import Fountain
from utils import xor, randChunkNums
import random


class DNAFountain(Fountain):
	def __init__(self, data=b'', chunk_size=35, index=1, seed=1, droplet_head_length=16, droplet_crc_length=16):
		#The length are DNA length in bps
		self.droplet_head_length = droplet_head_length
		#self.droplet_tail_length = droplet_tail_length
		self.droplet_crc_length = droplet_crc_length

		self.des = False
		self.sec_key = "hellhell"
		self.sec_key_int = 0
		super(DNAFountain, self).__init__(data, chunk_size, index, seed)


	def DNAdroplet(self):
		self.updateIndex()
		aDroplet = DNADroplet(b'', self.index, self.num_of_chunks, self.droplet_head_length,  self.seed, self.droplet_crc_length)
		aDroplet.sec_key = self.sec_key
		if self.des:
			aDroplet.des = True

			# print("DNAFountain sec_key: " + self.sec_key)

		aDroplet.data_len = self.chunk_size * 4
		aDroplet.degree = self.degrees[self.index % len(self.degrees)]
		chunk_nums = aDroplet.get_chunk_nums()
		data = None
		for num in chunk_nums:
			if data is None:
				data = self.chunk(num)
			else:
				data = xor(data, self.chunk(num))
		aDroplet.data = data
		if self.des:
			aDroplet.encry_data()
		return aDroplet

	def fix_bytes(self):
		blen = len(self.data)
		l_bytes = blen % self.chunk_size
		if l_bytes > 0:
			self.data = self.data + bytes(self.chunk_size - l_bytes)
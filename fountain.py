"""Shared chunk geometry and legacy soliton schedule for fountain encoding."""

from utils import get_degrees


class Fountain:
    def __init__(self, data=b'', chunk_size=35, index=1, seed=1):
        if type(chunk_size) is not int or chunk_size <= 0:
            raise ValueError("chunk size must be positive")
        if type(seed) is not int or seed <= 0:
            raise ValueError("fountain seed must be a positive integer")
        if type(index) is not int or index < 0:
            raise ValueError("initial droplet index must be a nonnegative integer")
        self.seed = seed
        self.index = index
        self.data = data
        self.original_size = len(data)
        self.chunk_size = chunk_size
        self.num_of_chunks = (len(data) + chunk_size - 1) // chunk_size
        self.ROBUST_FAILURE_PROBABILITY = 0.01
        self.c_value = 0.01
        self.degree_table_folds = 3
        self.VERBOSE = False
        self.gen_degrees()

    def chunk(self, num):
        start = self.chunk_size * num
        return self.data[start:start + self.chunk_size]

    def updateIndex(self):
        self.index += 1

    def updateSeed(self):
        self.seed += 1

    def gen_degrees(self):
        size = self.num_of_chunks * self.degree_table_folds
        if self.num_of_chunks < 2:
            self.degrees = [1] * size
        else:
            self.degrees = get_degrees(
                self.num_of_chunks, size, self.seed, self.ROBUST_FAILURE_PROBABILITY, self.c_value,
            )

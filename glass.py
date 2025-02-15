import sys
from utils import xor

sys.setrecursionlimit(100000)


class Glass:
    def __init__(self, num_chunks):
        self.entries = []
        self.droplets = []
        self.num_chunks = num_chunks
        self.chunks = [None] * num_chunks
        
    def addDroplet(self, d):
        self.droplets.append(d)
        entry = [d.get_chunk_nums(), d.data]
        self.entries.append(entry)
        #self.updateEntry(entry)

    def decode(self):
        if len(self.droplets) > 10:
            solved_blocks_count = 0
            iteration_solved_count = 0
            while iteration_solved_count > 0 or solved_blocks_count == 0:
                iteration_solved_count = 0
                for i, entry in enumerate(self.entries):
                    if len(entry[0]) == 1:
                        iteration_solved_count += 1
                        entry_index = next(iter(entry[0]))
                        # if self.chunks[entry_index] is not None:
                        #     assert self.chunks[entry_index] == entry[1], "Wrong droplet detected!"
                        self.chunks[entry_index] = entry[1]

                        self.entries.pop(i)

                        solved_blocks_count += 1
                        self.detach_entry(entry)
        else:
            print('Too few droplets!')

            # Try to recover more droplets
            # arr_in_arr
            # for i, entry_i in enumerate(self.entries):
            #     for j, entry_j in enumerate(self.entries):
            #         if arr_in_arr(entry_i[0], entry_j[0]) and len(entry_j[0]) - len(entry_i[0]) == 1:
            #             iteration_solved_count += 1
            #             entry_j[1] = xor(entry_j[1], entry_i[1])
            #             for e in entry_i[0]:
            #                 entry_j[0].remove(e)
            #             self.chunks[ next(iter(entry_j[0])) ] = entry_j[1]
            #             self.entries.pop(j)

    def detach_entry(self, entry):
        for other_entry in self.entries:
            if len(other_entry[0]) > 1 and next(iter(entry[0])) in other_entry[0]:
                other_entry[1] = xor(other_entry[1], entry[1])
                other_entry[0].remove(next(iter(entry[0])))


    def updateEntry(self, entry):
        for chunk_num in entry[0]:
            if self.chunks[chunk_num] is not None:
                entry[1] = xor(entry[1], self.chunks[chunk_num])
                entry[0].remove(chunk_num)
        if len(entry[0]) == 1:
            self.chunks[entry[0][0]] = entry[1]
            self.entries.remove(entry)
            for d in self.entries:
                if entry[0][0] in d[0]:
                    self.updateEntry(d)
                    
    def getString(self):
        return ''.join(x or ' _ ' for x in self.chunks)
        
    def isDone(self):
        return None not in self.chunks

    def chunksDone(self):
        count = 0
        for c in self.chunks:
            if c is not None:
                count+=1
        return count


    def toDNA(self):
        text = ''
        for drops in self.droplets:
            text = text + drops.toDNA() + "\n"
        return text


    def writeToFile(self, output_file):
        OUT = open(output_file, 'wb')
        i = 0
        while i < self.num_chunks:
            OUT.write(self.chunks[i])
            i = i + 1
        OUT.close()
        return True

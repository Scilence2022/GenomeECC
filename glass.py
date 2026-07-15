import logging

from utils import xor

logger = logging.getLogger(__name__)


class Glass:
    def __init__(self, num_chunks):
        self.entries = []
        self.droplets = []
        self.num_chunks = num_chunks
        self.chunks = [None] * num_chunks

    def addDroplet(self, d):
        self.droplets.append(d)
        # list(...) copies: the entry used to alias the droplet's own chunk_num_list, so
        # peeling entries during decode mutated the droplet's state as a side effect.
        chunk_nums = d.get_chunk_nums()
        entry = [list(chunk_nums) if chunk_nums else [], d.data]
        self.entries.append(entry)
        #self.updateEntry(entry)

    def decode(self):
        """Peel degree-1 droplets until no further progress is possible.

        Returns True if every chunk was recovered.
        """
        # Rewritten. The old loop condition (`while iteration_solved_count > 0 or
        # solved_blocks_count == 0`) spun forever when no degree-1 droplet was ever
        # present, it popped from self.entries while enumerating it (skipping entries),
        # and it refused to decode at all with <= 10 droplets.
        progress = True
        while progress:
            progress = False
            # Reverse order so pop(i) cannot shift an index we have yet to visit.
            for i in range(len(self.entries) - 1, -1, -1):
                entry = self.entries[i]
                if len(entry[0]) == 1:
                    entry_index = next(iter(entry[0]))
                    self.chunks[entry_index] = entry[1]
                    self.entries.pop(i)
                    self.detach_entry(entry)
                    progress = True
        if not self.isDone():
            logger.info("decode incomplete: %d/%d chunks recovered from %d droplets",
                        self.chunksDone(), self.num_chunks, len(self.droplets))
        return self.isDone()

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
        # chunks hold bytes, so ''.join() raised TypeError for any decoded glass.
        return b''.join(x if x else b' _ ' for x in self.chunks)
        
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
        if not self.isDone():
            raise ValueError(
                "cannot write output: %d of %d chunks are still unrecovered"
                % (self.num_chunks - self.chunksDone(), self.num_chunks))
        with open(output_file, 'wb') as OUT:
            for chunk in self.chunks:
                OUT.write(chunk)
        return True

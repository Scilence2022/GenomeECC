"""Incremental peeling decoder for verified fountain equations."""

from collections import deque
from hashlib import sha256
import logging

from ecc_format import write_atomic
from utils import xor


logger = logging.getLogger(__name__)


class Glass:
    def __init__(self, num_chunks, *, original_size=None, digest=None, chunk_size=None):
        if type(num_chunks) is not int or num_chunks < 0:
            raise ValueError("chunk count must be a nonnegative integer")
        if original_size is not None and (type(original_size) is not int or original_size < 0):
            raise ValueError("original payload size must be a nonnegative integer")
        if chunk_size is not None and (type(chunk_size) is not int or chunk_size <= 0):
            raise ValueError("chunk size must be positive")
        self.droplets = []
        self.num_chunks = num_chunks
        self.chunks = [None] * num_chunks
        self.original_size = original_size
        self.digest = digest
        self.chunk_size = chunk_size
        self._entries = {}
        self._by_chunk = [set() for _ in self.chunks]
        self._ready = deque()
        self._error = None

    def _check_error(self):
        if self._error:
            raise ValueError(self._error)

    def _check_residual(self, data):
        if any(data):
            self._error = "conflicting fountain equations"
            self._check_error()

    def addDroplet(self, d):
        self._check_error()
        indices = list(d.get_chunk_nums() or ())
        data = d.data
        if any(type(i) is not int or not 0 <= i < self.num_chunks for i in indices):
            raise ValueError("invalid fountain chunk index")
        if len(set(indices)) != len(indices):
            raise ValueError("duplicate fountain chunk index")
        if not isinstance(data, bytes) or (self.chunk_size is not None and len(data) != self.chunk_size):
            raise ValueError("invalid fountain payload width")
        pending = set()
        for i in indices:
            if self.chunks[i] is None:
                pending.add(i)
            else:
                data = xor(data, self.chunks[i])
        if not pending:
            self._check_residual(data)
            return
        entry_id = len(self.droplets)
        self.droplets.append(d)
        self._entries[entry_id] = (pending, data)
        for i in pending:
            self._by_chunk[i].add(entry_id)
        if len(pending) == 1:
            self._ready.append(entry_id)

    def decode(self):
        """Consume the work queue and return whether all chunks have been solved."""
        self._check_error()
        while self._ready:
            entry_id = self._ready.popleft()
            entry = self._entries.pop(entry_id, None)
            if entry is None:
                continue
            indices, data = entry
            i = next(iter(indices))
            self.chunks[i] = data
            self._by_chunk[i].discard(entry_id)
            for other_id in self._by_chunk[i]:
                remaining, payload = self._entries[other_id]
                remaining.remove(i)
                payload = xor(payload, data)
                if not remaining:
                    self._check_residual(payload)
                    del self._entries[other_id]
                    continue
                self._entries[other_id] = (remaining, payload)
                if len(remaining) == 1:
                    self._ready.append(other_id)
            self._by_chunk[i].clear()
        if not self.isDone():
            logger.info("decode incomplete: %d/%d chunks recovered", self.chunksDone(), self.num_chunks)
            return False
        if self.digest is not None or self.original_size is not None:
            try:
                self._complete_payload()
            except ValueError as exc:
                self._error = str(exc)
                raise
        return True

    def _complete_payload(self):
        self._check_error()
        if not self.isDone():
            raise ValueError(f"{self.num_chunks - self.chunksDone()} chunks are still unrecovered")
        data = b''.join(self.chunks)
        if self.original_size is not None:
            if len(data) < self.original_size or any(data[self.original_size:]):
                raise ValueError("invalid fountain padding or original size")
            data = data[:self.original_size]
        if self.digest is not None and sha256(data).hexdigest() != self.digest:
            raise ValueError("fountain payload checksum mismatch")
        return data

    def getString(self):
        self._check_error()
        if self.isDone():
            return self._complete_payload()
        return b''.join(chunk if chunk is not None else b' _ ' for chunk in self.chunks)

    def isDone(self):
        return not self._error and None not in self.chunks

    def chunksDone(self):
        return sum(chunk is not None for chunk in self.chunks)

    def toDNA(self):
        return ''.join(droplet.to_DNA() + '\n' for droplet in self.droplets)

    def writeToFile(self, output_file):
        write_atomic(output_file, self._complete_payload())
        return True

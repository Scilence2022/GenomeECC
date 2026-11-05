"""Manifest-based DNA fountain transport for binary payloads."""

from dataclasses import asdict, dataclass
from hashlib import sha256
import json

from DNAdroplet import DNADroplet
from DNAfountain import DNAFountain
from glass import Glass
from utils import get_degrees


@dataclass(frozen=True)
class FountainManifest:
    size: int
    chunk_size: int
    seed: int
    sha256: str
    version: str = "legacy-soliton-crc32-v1"

    def __post_init__(self):
        if self.version != "legacy-soliton-crc32-v1":
            raise ValueError("unsupported DNA fountain format version")
        if type(self.size) is not int or self.size < 0:
            raise ValueError("fountain size must be a nonnegative integer")
        if type(self.chunk_size) is not int or self.chunk_size <= 0:
            raise ValueError("fountain chunk size must be positive")
        if type(self.seed) is not int or self.seed <= 0:
            raise ValueError("fountain seed must be a positive integer")
        if type(self.sha256) is not str or len(self.sha256) != 64 or any(
            char not in "0123456789abcdef" for char in self.sha256
        ):
            raise ValueError("invalid fountain SHA-256")

    @property
    def chunks(self):
        return (self.size + self.chunk_size - 1) // self.chunk_size

    def to_bytes(self):
        return json.dumps(asdict(self), sort_keys=True).encode('ascii')

    @classmethod
    def from_bytes(cls, data):
        try:
            fields = json.loads(data)
            if type(fields) is not dict or fields.keys() != {
                "size", "chunk_size", "seed", "sha256", "version",
            }:
                raise ValueError("invalid fountain manifest fields")
            return cls(**fields)
        except (UnicodeError, RecursionError) as exc:
            raise ValueError("invalid fountain manifest") from exc


def encode_dna(data, count, *, chunk_size=35, seed=1, strand_filter=None, max_attempts=None):
    manifest = FountainManifest(len(data), chunk_size, seed, sha256(data).hexdigest())
    if type(count) is not int or count < manifest.chunks:
        raise ValueError("droplet count must be at least the number of chunks")
    if not manifest.size:
        return manifest, []
    budget = count * 10 if max_attempts is None else max_attempts
    if type(budget) is not int or budget < count:
        raise ValueError("attempt limit must be at least the requested droplet count")
    fountain = DNAFountain(data, chunk_size=chunk_size, seed=seed)
    strands = []
    for _ in range(budget):
        strand = fountain.DNAdroplet().to_DNA()
        if strand_filter is None or strand_filter(strand):
            strands.append(strand)
            if len(strands) == count:
                return manifest, strands
    raise ValueError(f"droplet attempt budget exhausted: accepted {len(strands)}/{count}")


def decode_dna(manifest, strands):
    chunks = manifest.chunks
    decoder = Glass(
        chunks, original_size=manifest.size, digest=manifest.sha256, chunk_size=manifest.chunk_size,
    )
    degrees = get_degrees(chunks, chunks * 3, manifest.seed, .01, .01) if chunks >= 2 else [1] * 3
    if chunks:
        for strand in strands:
            droplet = DNADroplet(bytes(manifest.chunk_size), num_of_chunks=chunks)
            if not droplet.set_droplet_from_DNA(strand):
                continue
            droplet.degree = degrees[droplet.index % len(degrees)]
            decoder.addDroplet(droplet)
        decoder.decode()
    if not decoder.isDone():
        raise ValueError(f"{chunks - decoder.chunksDone()} fountain chunks are still unrecovered")
    return decoder.getString()

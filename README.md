# GenomeECC

GenomeECC creates separate error-correction archives for DNA sequences. It packs
A/C/G/T into two bits per base, stores anchors and per-block checksums, and applies
Reed–Solomon protection across columns of sequence blocks.

## Installation

```sh
python -m pip install -e .
python -m pip install -e '.[dev]'
```

Core encoding, verification, and recovery from known blocks require Python 3.9
or newer. Optional read assembly uses khmer. A tested configuration uses the
Bioconda khmer 3.0.0a3 build with Python 3.11.

```sh
conda create -n genomeecc-assembly --override-channels -c conda-forge -c bioconda \
  --strict-channel-priority python=3.11 khmer=3.0.0a3
conda activate genomeecc-assembly
python -m pip install -e '.[dev]'
```

[Bioconda builds by platform](https://anaconda.org/channels/bioconda/packages/khmer/files)

## Encode and verify

```sh
genome-ecc encode -g genomes/cgl-ssbl-from_nxm.fna -o genome.ecc \
  --block-len 320 --chunk-size 10 --redundancy 0.04
genome-ecc info -e genome.ecc
genome-ecc verify -e genome.ecc -g genomes/cgl-ssbl-from_nxm.fna
genome-ecc verify -e genome.ecc -g genomes/cgl-ssbl-from_nxm.fna --exact
```

The default v2 archive records FASTA IDs, titles, topology, per-record SHA-256,
and the boundaries of independent RS groups. `info` validates the file layout
and the v2 archive checksum. `verify` compares the reference sequences with the
stored digests and uses the archived parameters. With `--exact`, verification
also regenerates the stored group bytes. Legacy records are matched by order,
and historical padding differences can cause an exact comparison to fail after
sequence verification succeeds.

FASTA input accepts upper- and lowercase A/C/G/T. Output normalizes bases to
uppercase. Duplicate IDs, empty records, and other nucleotide symbols produce an
error with input context. Put `topology=circular` or `topology=linear` in the
record description to specify topology explicitly.

## Recover

Given positioned blocks in a tab-separated file, run the following command.

```sh
genome-ecc recover -e genome.ecc --blocks candidates.tsv -o recovered.fna
```

Each TSV row contains the record ID, a zero-based block index within that record,
and the DNA sequence. Use zero-based record numbers for legacy files. Leave out
missing blocks or write `-` in the sequence column. Supply full block lengths,
with the actual shorter length accepted for the final block. CRC failures count
as erasures. The command replaces the destination after all records pass their
sequence digest checks and reports missing and CRC-rejected blocks.

For FASTQ or FASTA reads, install the optional k-mer backend and use `--reads`.
The assembler searches from stored anchors and passes unambiguous candidates to
the same recovery code. Legacy archives use linear topology for assembly.

```sh
genome-ecc recover -e genome.ecc --reads reads_R1.fq reads_R2.fq \
  --kmer-size 31 --min-kmer-coverage 5 --max-paths 1000 -o recovered.fna
```

The default k-mer length is 31 or the shortest anchor width in bases, whichever
is smaller. Insufficient coverage,
repeats, and indels can leave blocks unresolved. The archive must have enough
parity per group for the resulting errors and erasures.

## Choose parameters

| Option | Default | Meaning |
|---|---:|---|
| `--block-len` | 320 | Bases per block. Use a multiple of 4 in 4–1020. A multiple of 8 avoids column padding. |
| `--chunk-size` | 10 | Blocks per anchor and CRC chunk, in 1–63. |
| `--anchor-bytes` | 8 | Bytes per anchor, in 0–15 and within chunk capacity. Zero is useful for positioned block recovery. |
| `--ecc-mode` | `crc8` | Per-block CRC mode. `crc16` uses two bytes and reduces accidental matches. |
| `--redundancy` | 0.04 | Parity rate in 0–1. Zero provides checksum verification without recovery capacity. |
| `--group-blocks` | 4096 | v2 group cap, rounded down to a whole number of chunks. |

Every group contains one GF(2¹⁶) codeword per column. A data block occupies one
symbol in each column. A group with `B` blocks and `R = ceil(rate × B)` parity
symbols can correct `E` unknown symbol errors and `S` known erasures per column
when `2E + S ≤ R`. The budget is separate for every group. A long contiguous
deletion can exhaust one group's parity even if other groups have unused parity.
A legacy single-codeword record requires `B + R ≤ 65,535`.

On the repository's 3.28 Mbp sample, the default v2 configuration produced a
51,631-byte archive in 7.4 seconds on one Apple M3 Max. The same parameters in
the legacy genome layout produced 51,378 bytes in 19.7 seconds. Recovery from
positioned blocks with 75 erased blocks took 21.2 seconds. These are single-run
wall-clock measurements with Python 3.9.6 and reedsolo 1.7.0.

## Legacy compatibility

The CLI detects legacy genome containers and naked chromosome ECC files by
checking the entire layout. Use `--format legacy-genome` or
`--format legacy-chromosome` to choose an explicit legacy format. The
`ChromosomeECC` and `GenomeECC` compatibility APIs also write legacy layouts.

```sh
genome-ecc info -e genomes/cgl-ssbl-from_nxm.fna.ecc
genome-ecc migrate -e genomes/cgl-ssbl-from_nxm.fna.ecc \
  -g genomes/cgl-ssbl-from_nxm.fna -o genome-v2.ecc
```

Migration checks the supplied reference against the historical digest and creates
a new v2 archive with the selected encoding parameters. Historical files whose
lengths disagree with their geometry, including old encodings that left a column
unprotected, are rejected. Such archives need a trusted source sequence for
re-encoding. The [binary format description](docs/ecc-format.md) specifies the
fields and integrity checks.

## DNA Fountain

Use the manifest API for a byte-exact round trip from DNA strand strings.

```python
from fountain_codec import FountainManifest, decode_dna, encode_dna

payload = b"example\x00"
manifest, strands = encode_dna(payload, count=20, chunk_size=4, seed=11)
manifest = FountainManifest.from_bytes(manifest.to_bytes())
restored = decode_dna(manifest, strands)
assert restored == payload
```

The manifest records the source length, SHA-256, chunk size, seed, and schedule
version. Decoding validates each strand's index and CRC32 and verifies the exact
source bytes. The lower-level Glass decoder also supports incremental batches.
Use a sufficient strand surplus
for the intended loss model. `strand_filter` and `max_attempts` allow bounded
screening for experimental sequence constraints.

The versioned schedule preserves the historical nonstandard soliton weighting.
It gives degree 1 additional mass, omits degree `K` for `K ≥ 2`, and places
the robust spike one degree below the original formula. Existing seed vectors
are covered by regression tests. [Original LT Codes paper](https://www.mit.edu/~6.454/www_fall_2003/khisti/luby02lit.pdf)

Archived experimental scripts such as `Encode.py`, `decode.py`, and
`genGenomeECC*.py` contain machine-specific inputs and older APIs. Use the
supported CLI and manifest API for new runs. The opt-in DES code is retained for
legacy droplet compatibility. It uses a key-derived permutation and the key as
its CBC IV. Apply authenticated encryption to the source payload before fountain
encoding if confidentiality or authenticity is required.

## Tests

```sh
python -m pytest
```

The suite covers the core code, both legacy layouts, v2 groups, on-disk recovery,
historical fixtures, malformed inputs, forward/reverse assembly, and real DNA
fountain strings. The khmer integration test runs when khmer is installed. CI
builds the distribution and runs the installed CLI outside the source tree.

## License

GPL v3. See [LICENSE](LICENSE).

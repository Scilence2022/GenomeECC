# GenomeECC

Error-correcting codes (ECC) for **whole-genome DNA data preservation**.

A genome stored in a living cell, or synthesised as DNA, accumulates damage: bases
mutate, stretches drop out, reads fail. GenomeECC generates a small, separate ECC file
that lets a damaged genome be reconstructed exactly. The ECC file for a 3.3 Mbp genome
is ~50 KB (about 6% of the packed sequence) and can restore whole erased blocks.

There is also a **DNA fountain code** subsystem (Luby-transform droplets over DNA) for
general-purpose DNA data storage, plus a De Bruijn graph assembler for recovering
droplets from sequencing reads.

---

## Install

```bash
pip install -e .                  # ECC pipeline + CLI
pip install -e '.[dev]'           # + pytest
pip install -e '.[assembly]'      # + khmer, for deBruijnGraph (Python <= 3.10 only)
```

`khmer` has no wheels for Python >= 3.11 and fails to build there. It is required
**only** by `deBruijnGraph.py`; the ECC pipeline does not need it.

## Quick start

```bash
# Generate ECC codes for a genome
genome-ecc encode -g genomes/cgl-ssbl-from_nxm.fna -o genome.ecc \
    --block-len 320 --chunk-size 10 --redundancy 0.04

# Describe an existing ECC file
genome-ecc info -e genome.ecc

# Check an ECC file still matches its genome
genome-ecc verify -g genomes/cgl-ssbl-from_nxm.fna -e genome.ecc \
    --block-len 320 --chunk-size 10 --redundancy 0.04
```

As a library:

```python
from chromosomeECC import ChromosomeECC

c = ChromosomeECC(block_len=320, chunk_size=10, redundant_block_rate=0.04)
c.set_chr_seq(sequence)             # sequence is an ACGT string
ecc_bytes = c.get_all_ecc_bytes()   # store this alongside the genome
```

## How it works

```
genome (ACGT)
  -> 2-bit packing            4 bases per byte (A=00 C=01 G=10 T=11)
  -> blocks                   block_len bases each (block_size = block_len/4 bytes)
  -> chunks                   chunk_size blocks each, each with an anchor + per-block CRCs
  -> columns                  byte i of every block forms column i
  -> Reed-Solomon per column  in GF(2^16), 2-byte symbols, one symbol per block
```

The column transposition is the important idea: because each column holds exactly **one
symbol per block**, losing an entire block costs every column exactly one symbol. A
block erasure therefore becomes a single-symbol erasure in each column, which
Reed-Solomon corrects independently. With `redundant_block_rate = r`, up to
`ceil(r * block_num)` whole blocks can be erased and rebuilt exactly.

The ECC file contains: genome length, MD5, block/chunk geometry, per-chunk anchors,
per-block CRCs, and the RS parity for each column.

## Choosing parameters

| Parameter | Meaning | Notes |
|---|---|---|
| `--block-len` | bases per block | 4..1020 (10-bit field). A multiple of 8 keeps `block_size` even and avoids column padding. |
| `--chunk-size` | blocks per chunk | 1..63 (6-bit field). |
| `--redundancy` | RS redundant block rate | Directly sets how many whole blocks can be erased and recovered. |
| `--ecc-mode` | per-block checksum | `crc8` (1 byte) or `crc16` (2 bytes). |

Values outside those ranges now raise `ValueError` rather than silently corrupting the
header.

## Testing

```bash
pytest
```

97 tests cover the encode/decode round trip, RS erasure recovery at and beyond the
correction limit, header packing, FASTA parsing, the codec, and the fountain decoder.

## Performance

On a real 3.28 Mbp genome (`genomes/cgl-ssbl-from_nxm.fna`, `block_len=320`):

| Stage | Before | After |
|---|---|---|
| `set_chr_seq` (DNA->bytes, CRCs, anchors) | 12.75 s | **0.14 s** |
| RS encode (third-party `reedsolo`) | 19.84 s | 19.87 s |
| **Total** | **32.60 s** | **20.03 s** |

`DNAToBytes` previously rebuilt a 256-entry lookup dict on every 4-base conversion; it
is now a vectorised 2-bit pack (~2400x faster in isolation). The remaining runtime is
almost entirely inside `reedsolo`'s pure-Python `rs_encode_msg`. If RS encoding
dominates your workload, that library is the place to optimise.

## Compatibility with previously generated .ECC files

For a configuration the old code handled correctly (e.g. `block_len=320`, explicit
`set_block_ecc_mode('crc8')`), the ECC output is unchanged except for **the final
chunk's unused CRC padding**. A short final chunk has fewer real blocks than
`chunk_size`; the old code filled the leftover slots by cyclically repeating earlier
CRCs (a side effect of `np.resize`, which tiles rather than pads), so those slots held
checksums for blocks that do not exist. They are now zero-filled. Reed-Solomon parity is
byte-for-byte identical.

Output **does** legitimately change where the old behaviour was broken — most
importantly when `block_len` gives an odd `block_size` (the defaults 500 and 300 did),
because the last byte of every block previously belonged to no column and was neither
protected nor restored.

## Known deviations from the published algorithms

`utils.py` implements the Luby soliton distributions with **deliberate deviations**,
preserved as-is so previously tuned experimental results stay reproducible. They are
documented here rather than silently "corrected":

- `ideal_distribution` uses `rho(1) = 10/K` instead of Luby's `1/K`, putting roughly 5x
  more mass on degree 1 after normalisation, and `range(2, K)` means degree `K` is
  never emitted.
- `robust_distribution` places the tau spike at degree `M-1` rather than `M`, and tau
  covers degrees `1..M-2` rather than `1..M-1`.

Both still return correctly normalised distributions. If you want textbook Luby
behaviour, change these two functions deliberately and re-run your experiments.

## Security notes

**Do not rely on this codebase for confidentiality.** The optional droplet encryption
(`des=True`) is broken by design and kept only for compatibility:

- It uses **DES**, which has a 56-bit key and has been considered broken since 1998.
- `des_en`/`des_de` set the CBC **IV equal to the key**, making encryption
  deterministic.
- The byte permutation in `DNAdroplet` derives its seed from `index + key`, while
  `index` is transmitted in plaintext and the CRC is computed over the ciphertext. That
  combination is an unkeyed oracle: an attacker holding only a public droplet can
  confirm a guessed key without any known plaintext. A working key-recovery
  demonstration was produced during a security review of this repository.

If you need real confidentiality, encrypt the payload with a modern AEAD (for example
`cryptography`'s `AESGCM`) before it reaches this library, and never derive
non-secret-looking state from the key.

## Repository layout

| File | Purpose |
|---|---|
| `chromosomeECC.py` | Core ECC for one chromosome (blocks, columns, RS, CRCs, anchors) |
| `genomeECC.py` | Multi-chromosome wrapper, FASTA reading, ECC file packing |
| `chromosomeChunk.py` | A chunk of blocks and its checksums |
| `genome_ecc_cli.py` | `genome-ecc` command-line interface |
| `utils.py` | DNA<->bytes codec, soliton distributions, GC/homopolymer checks |
| `test_utils.py` | Helper/experiment library (**not** a test module, despite the name) |
| `deBruijnGraph.py` | De Bruijn graph assembly; the only user of the optional `khmer` |
| `DNAdroplet.py`, `DNAfountain.py`, `droplet.py`, `fountain.py`, `glass.py` | Fountain-code subsystem |
| `genGenomeECC*.py` | Original experiment scripts, kept as a record. They contain absolute paths to the author's machine and are not runnable elsewhere -- use `genome-ecc` instead. |
| `tests/` | pytest suite |

## License

GPL v3. See [LICENSE](LICENSE).

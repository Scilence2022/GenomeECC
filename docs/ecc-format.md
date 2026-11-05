# GenomeECC archive layouts

Integers use unsigned big-endian byte order. GenomeECC packs four bases in each
byte, using A=0, C=1, G=2, and T=3 with the first base in the high bits.

## Version 2 envelope

| Field | Size | Meaning |
|---|---:|---|
| Magic | 8 bytes | `47 45 43 43 02 0d 0a 1a` |
| Metadata length | 4 bytes | Length of the UTF-8 JSON section, up to 16 MiB |
| Metadata | Variable | Nonempty JSON array of records, in FASTA order |
| Group payloads | Variable | Concatenated chromosome payloads in metadata order |
| Archive digest | 32 bytes | SHA-256 of every preceding byte |

Each JSON record contains `id`, `title` (the FASTA title without the leading
`>`), `topology` (`linear`, `circular`, or `unknown`), `sha256` (lowercase hex
of the uppercase ASCII bases), and `groups` (the byte length of each payload).
IDs are unique. The ordered group lengths delimit payloads and must consume
the archive exactly. Each group covers the next contiguous range of source
blocks. Nonterminal groups end at a chunk boundary. Groups in a record use
the same block size, chunk size, anchor width, and CRC mode. Their parity
counts can differ.

The format stores both a full-record sequence SHA-256 and the historical group
MD5s. A reader checks the envelope digest before using the metadata. An external
signature or authenticated copy is needed for authenticity against an active
attacker.

## Chromosome payload

A v2 group, a standalone legacy chromosome, and each record inside a legacy
genome container use the same 25-byte payload header.

| Field | Size | Meaning |
|---|---:|---|
| Bases | 4 bytes | True length before padding |
| Group digest | 16 bytes | MD5 of the 2-bit sequence padded to a multiple of four bases with A |
| Geometry | 2 bytes | High 10 bits are bases per block, low 6 bits are blocks per chunk |
| Widths | 1 byte | High nibble is anchor width, low nibble is block CRC width, in bytes |
| Parity count | 2 bytes | Number of RS parity symbols per column |

The header is followed by one fixed-width metadata entry per chunk. Each entry
holds its anchor, padded with zero bytes, then one CRC per slot. The final
short block is padded with A for computing its CRC. New archives zero-fill
unused CRC slots. The CRC algorithms are CRC-8 with polynomial `0x07`, and
CRC-16/XMODEM with polynomial `0x1021`.

The remaining section contains the RS parity symbols of each column in order.
Each block is padded to an even number of bytes. The first two bytes from every
block form column zero, the next two bytes form column one, and so on. One
column uses a single Reed–Solomon codeword over GF(2¹⁶) with polynomial
`0x1002d`, generator 2, and first consecutive root 0. Symbols use two bytes.
Every codeword must contain at most 65,535 symbols including parity.

For `B` blocks, `M` blocks per chunk, `A` anchor bytes, `C` CRC bytes, `W`
packed bytes per block, and `R` parity symbols per column, the payload byte
length is `25 + ceil(B/M) × (A + M × C) + 2 × ceil(W/2) × R`.

## Legacy genome container

A one-byte record count is followed by a four-byte payload length and a
chromosome payload for each record. Record IDs and topology are external to
this layout. Standalone legacy chromosome files contain only a chromosome
payload. Readers check each candidate layout fully before selecting it.

The reader accepts the historical MD5 convention that appended an extra byte
of zero to already aligned DNA, and skips unused CRC slots. It requires the
declared widths and parity to cover the full layout. Migration takes a matching
reference FASTA and writes a fresh v2 archive.

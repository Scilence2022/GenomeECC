#!/usr/bin/env python3
"""Command-line interface for generating and inspecting genome ECC codes.

The genGenomeECC*.py scripts are kept as a record of specific experiments and have
absolute paths baked in. This module does the same job with explicit arguments so it
runs on any machine:

    genome-ecc encode -g genome.fna -o genome.ecc --block-len 320 --chunk-size 10
    genome-ecc info   -e genome.ecc
    genome-ecc verify -g genome.fna -e genome.ecc
"""
import argparse
import logging
import sys

from chromosomeECC import ChromosomeECC
from genomeECC import GenomeECC

logger = logging.getLogger("genome_ecc")


def _build(args):
    g = GenomeECC(block_len=args.block_len, chunk_size=args.chunk_size,
                  redundant_block_rate=args.redundancy)
    g.read_genome(args.genome)
    if not g.chr_num:
        raise SystemExit("error: no FASTA records found in %s" % args.genome)
    for head in g.genome_seqs:
        chr_ecc = g.genChromesomeECC(head, g.genome_seqs[head])
        chr_ecc.set_block_ecc_mode(args.ecc_mode)
        chr_ecc.chunk_anchor_byte_num = args.anchor_bytes
        chr_ecc.set_chr_seq(g.genome_seqs[head])
        g.chr_ECC_objs[head] = chr_ecc
    return g


def _split_chr_blobs(ecc_bytes):
    """Yield each per-chromosome ECC blob from a genome-level ECC file."""
    chr_num = int.from_bytes(ecc_bytes[:1], 'big', signed=False)
    pos = 1
    for _ in range(chr_num):
        blob_len = int.from_bytes(ecc_bytes[pos:pos + 4], 'big', signed=False)
        yield ecc_bytes[pos + 4:pos + 4 + blob_len]
        pos += 4 + blob_len


def cmd_encode(args):
    g = _build(args)
    blob = g.get_all_genome_ecc()
    with open(args.output, 'wb') as fh:
        fh.write(blob)

    genome_bases = sum(len(s) for s in g.genome_seqs.values())
    # 4 bases pack into 1 byte, so this is ECC size relative to the packed genome.
    packed = genome_bases / 4 or 1
    print("chromosomes : %d" % g.chr_num)
    print("total bases : %d" % genome_bases)
    print("ECC bytes   : %d" % len(blob))
    print("overhead    : %.2f%% of packed genome" % (100.0 * len(blob) / packed))
    print("written to  : %s" % args.output)
    return 0


def cmd_info(args):
    with open(args.ecc, 'rb') as fh:
        blob = fh.read()
    blobs = list(_split_chr_blobs(blob))
    print("ECC file    : %s" % args.ecc)
    print("total bytes : %d" % len(blob))
    print("chromosomes : %d" % len(blobs))
    probe = ChromosomeECC()
    for i, chr_blob in enumerate(blobs):
        h = probe.unpack_ecc_bytes(chr_blob)
        print("  [%d] length=%d block_len=%d chunk_size=%d blocks=%d chunks=%d "
              "anchor=%dB crc=%dB redundancy=%d md5=%s"
              % (i, h['genome_length'], h['block_len'], h['chunk_size'],
                 h['block_num'], h['chunk_num'], h['chunk_anchor_byte_num'],
                 h['block_crc_byte_num'], h['redundant_block_num'],
                 h['genome_md5'].hex()[:16]))
    return 0


def cmd_verify(args):
    """Re-derive ECC from the FASTA and check it matches the stored file."""
    with open(args.ecc, 'rb') as fh:
        stored = fh.read()
    g = _build(args)
    fresh = g.get_all_genome_ecc()
    if fresh == stored:
        print("OK: regenerated ECC is byte-identical to %s" % args.ecc)
        return 0
    print("MISMATCH: regenerated ECC differs from %s" % args.ecc, file=sys.stderr)
    print("  stored=%d bytes regenerated=%d bytes" % (len(stored), len(fresh)),
          file=sys.stderr)
    print("  (differing parameters produce a different file -- check --block-len, "
          "--chunk-size, --redundancy, --ecc-mode, --anchor-bytes)", file=sys.stderr)
    return 1


def _add_encode_args(p):
    p.add_argument('-g', '--genome', required=True, help="input genome FASTA")
    p.add_argument('--block-len', type=int, default=320,
                   help="block length in bases, 4..1020, multiple of 8 recommended "
                        "(default: 320)")
    p.add_argument('--chunk-size', type=int, default=10,
                   help="blocks per chunk, 1..63 (default: 10)")
    p.add_argument('--redundancy', type=float, default=0.04,
                   help="Reed-Solomon redundant block rate (default: 0.04)")
    p.add_argument('--ecc-mode', choices=['crc8', 'crc16'], default='crc8',
                   help="per-block checksum width (default: crc8)")
    p.add_argument('--anchor-bytes', type=int, default=8,
                   help="chunk anchor length in bytes (default: 8)")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='genome-ecc',
        description="Generate and inspect error-correcting codes for genome preservation.")
    parser.add_argument('-v', '--verbose', action='store_true', help="debug logging")
    sub = parser.add_subparsers(dest='command', required=True)

    p_enc = sub.add_parser('encode', help="generate ECC codes for a genome")
    _add_encode_args(p_enc)
    p_enc.add_argument('-o', '--output', required=True, help="output .ecc file")
    p_enc.set_defaults(func=cmd_encode)

    p_info = sub.add_parser('info', help="describe an existing .ecc file")
    p_info.add_argument('-e', '--ecc', required=True, help="input .ecc file")
    p_info.set_defaults(func=cmd_info)

    p_ver = sub.add_parser('verify', help="check an .ecc file against its genome")
    _add_encode_args(p_ver)
    p_ver.add_argument('-e', '--ecc', required=True, help="the .ecc file to check")
    p_ver.set_defaults(func=cmd_verify)

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")
    try:
        return args.func(args)
    except (OSError, ValueError) as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

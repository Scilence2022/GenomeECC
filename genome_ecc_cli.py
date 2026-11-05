#!/usr/bin/env python3
"""Encode, inspect, verify, and recover GenomeECC archives."""

import argparse
import logging
from pathlib import Path
import sys

from ecc_archive import FORMATS, encode_archive, parse_archive, recover_blocks, verify_records
from ecc_format import CRC_MODES, ceil_div, write_atomic
from fasta import read_fasta


def _load_archive(args):
    blob = Path(args.ecc).read_bytes()
    return parse_archive(blob, args.format), len(blob)


def _check_output(output, *inputs):
    resolved = Path(output).resolve()
    if any(resolved == Path(source).resolve() for source in inputs):
        raise ValueError("output must differ from the input files")


def _encode(args, records, format):
    blob = encode_archive(
        records, block_len=args.block_len, chunk_size=args.chunk_size,
        rate=args.redundancy, mode=args.ecc_mode, anchor=args.anchor_bytes,
        format=format, group_blocks=args.group_blocks,
    )
    write_atomic(args.output, blob)
    return len(blob)


def cmd_encode(args):
    _check_output(args.output, args.genome)
    bases = packed = count = 0

    def records():
        nonlocal bases, packed, count
        for record in read_fasta(args.genome):
            bases += len(record.sequence)
            packed += ceil_div(len(record.sequence), 4)
            count += 1
            yield record

    size = _encode(args, records(), args.format)
    print(f"format      : {args.format}")
    print(f"chromosomes : {count}")
    print(f"total bases : {bases}")
    print(f"ECC bytes   : {size}")
    print(f"overhead    : {size / packed:.2%} of packed genome")
    print(f"written to  : {args.output}")
    return 0


def cmd_info(args):
    archive, size = _load_archive(args)
    print(f"ECC file    : {args.ecc}")
    print(f"format      : {archive.format}")
    print(f"total bytes : {size}")
    print(f"chromosomes : {len(archive.records)}")
    for i, record in enumerate(archive.records):
        header = record.groups[0].header
        parity = ','.join(str(group.header['redundant_block_num']) for group in record.groups)
        digest = f"sha256={record.digest}" if record.digest else f"md5={header['genome_md5'].hex()}"
        print(
            f"  [{i}] id={record.name} topology={record.topology} length={record.length} "
            f"block_len={header['block_len']} chunk_size={header['chunk_size']} "
            f"blocks={record.blocks} groups={len(record.groups)} "
            f"anchor={header['chunk_anchor_byte_num']}B crc={header['block_crc_byte_num']}B "
            f"parity_per_group={parity} {digest}"
        )
    return 0


def cmd_verify(args):
    archive, _ = _load_archive(args)
    verify_records(
        archive, read_fasta(args.genome), exact=args.exact, block_len=args.block_len,
        chunk_size=args.chunk_size, rate=args.redundancy, mode=args.ecc_mode, anchor=args.anchor_bytes,
    )
    result = "regenerated ECC matches" if args.exact else "sequence checksums match"
    print(f"OK: {result} for {len(archive.records)} record(s)")
    return 0


def _read_blocks(path):
    candidates = {}
    with open(path, encoding="utf-8") as file:
        for number, line in enumerate(file, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split('\t')
            if len(fields) != 3 or not all(fields) or not fields[1].isascii() or not fields[1].isdigit():
                raise ValueError(f"{path}:{number}: expected record<TAB>block_index<TAB>sequence")
            record, index, sequence = fields
            index = int(index)
            blocks = candidates.setdefault(record, {})
            if index in blocks:
                raise ValueError(f"{path}:{number}: duplicate block {record!r}/{index}")
            blocks[index] = None if sequence == "-" else sequence
    return candidates


def cmd_recover(args):
    sources = [args.blocks] if args.blocks else args.reads
    _check_output(args.output, args.ecc, *sources)
    archive, _ = _load_archive(args)
    if args.blocks:
        candidates = _read_blocks(args.blocks)
    else:
        from deBruijnGraph import DeBruijnGraph, assemble_candidates

        if args.min_kmer_coverage < 1 or args.max_paths < 1:
            raise ValueError("k-mer coverage and maximum paths must be positive")
        anchors = [
            4 * record.groups[0].header["chunk_anchor_byte_num"] for record in archive.records
            if record.length > 4 * record.groups[0].header["chunk_anchor_byte_num"]
        ]
        anchor_limit = min(anchors, default=32)
        kmer_size = min(31, anchor_limit) if args.kmer_size is None else args.kmer_size
        if not 1 <= kmer_size <= anchor_limit:
            raise ValueError(
                f"k-mer size must fit the available anchors (maximum {anchor_limit}); "
                "use --blocks for archives without anchors"
            )
        graph = DeBruijnGraph(kmer_size)
        graph.min_kmer_cov = args.min_kmer_coverage
        graph.max_path_num = args.max_paths
        for reads in args.reads:
            graph.count_file(reads)
        candidates = assemble_candidates(archive, graph)
    results = recover_blocks(archive, candidates)
    def fasta_lines():
        for record in results:
            yield f">{record.title}\n".encode("utf-8")
            for i in range(0, len(record.sequence), 80):
                yield record.sequence[i:i + 80].encode("ascii") + b'\n'

    write_atomic(args.output, fasta_lines())
    for record in results:
        print(
            f"{record.name}: {len(record.sequence)} bases verified, "
            f"{record.missing_blocks} missing blocks, {record.rejected_blocks} CRC-rejected blocks, "
            f"{record.corrected_symbols} corrected RS symbols"
        )
    return 0


def cmd_migrate(args):
    _check_output(args.output, args.genome, args.ecc)
    archive, _ = _load_archive(args)
    records = list(read_fasta(args.genome))
    verify_records(archive, records)
    size = _encode(args, records, "v2")
    print(f"OK: wrote verified v2 archive ({size} bytes) to {args.output}")
    return 0


def _add_geometry(parser, *, encoding):
    defaults = (320, 10, .04, "crc8", 8) if encoding else (None,) * 5
    parser.add_argument('--block-len', type=int, default=defaults[0],
                        help="bases per block, a multiple of 4 in 4..1020")
    parser.add_argument('--chunk-size', type=int, default=defaults[1],
                        help="blocks per chunk, 1..63")
    parser.add_argument('--redundancy', type=float, default=defaults[2],
                        help="Reed–Solomon parity rate in 0..1")
    parser.add_argument('--ecc-mode', choices=tuple(CRC_MODES.values()), default=defaults[3])
    parser.add_argument('--anchor-bytes', type=int, default=defaults[4], help="anchor width in 0..15")
    if encoding:
        parser.add_argument('--group-blocks', type=int, help="v2 group cap, rounded down to full chunks")


def _add_archive_input(parser):
    parser.add_argument('-e', '--ecc', required=True, help="input ECC archive")
    parser.add_argument('--format', choices=('auto', *FORMATS), default='auto')


def main(argv=None):
    parser = argparse.ArgumentParser(prog="genome-ecc", description=__doc__)
    parser.add_argument('-v', '--verbose', action='store_true', help="debug logging")
    sub = parser.add_subparsers(dest='command', required=True)

    encode = sub.add_parser('encode', help="generate an ECC archive")
    encode.add_argument('-g', '--genome', required=True, help="input genome FASTA")
    encode.add_argument('-o', '--output', required=True, help="output ECC archive")
    encode.add_argument('--format', choices=FORMATS, default='v2')
    _add_geometry(encode, encoding=True)
    encode.set_defaults(func=cmd_encode)

    info = sub.add_parser('info', help="validate structure and describe an ECC archive")
    _add_archive_input(info)
    info.set_defaults(func=cmd_info)

    verify = sub.add_parser('verify', help="check sequence checksums against an archive")
    _add_archive_input(verify)
    verify.add_argument('-g', '--genome', required=True, help="reference genome FASTA")
    verify.add_argument('--exact', action='store_true', help="also compare regenerated group bytes")
    _add_geometry(verify, encoding=False)
    verify.set_defaults(func=cmd_verify)

    recover = sub.add_parser('recover', help="recover and verify a FASTA")
    _add_archive_input(recover)
    inputs = recover.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--blocks', help="TSV with record ID, zero-based block index, and bases")
    inputs.add_argument('--reads', nargs='+', help="FASTA/FASTQ inputs for optional k-mer assembly")
    recover.add_argument('--kmer-size', type=int, help="k-mer size, capped by the shortest anchor")
    recover.add_argument('--min-kmer-coverage', type=int, default=5)
    recover.add_argument('--max-paths', type=int, default=1000)
    recover.add_argument('-o', '--output', required=True, help="verified output FASTA")
    recover.set_defaults(func=cmd_recover)

    migrate = sub.add_parser('migrate', help="verify a reference and encode a new v2 archive")
    _add_archive_input(migrate)
    migrate.add_argument('-g', '--genome', required=True, help="matching reference genome FASTA")
    migrate.add_argument('-o', '--output', required=True, help="new v2 archive")
    _add_geometry(migrate, encoding=True)
    migrate.set_defaults(func=cmd_migrate)

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")
    try:
        return args.func(args)
    except (OSError, ValueError, ImportError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

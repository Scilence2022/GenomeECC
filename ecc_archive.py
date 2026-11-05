"""Versioned ECC archives, historical compatibility, and verified block recovery."""

from dataclasses import dataclass
from hashlib import md5, sha256
import json
import math

from reedsolo import ReedSolomonError

from chromosomeChunk import block_crc
from chromosomeECC import ChromosomeECC
from ecc_format import (
    CRC_MODES, MAX_CODEWORD, EccFormatError, ceil_div, unpack_chromosome,
    unpack_legacy_genome, validate_codeword, validate_geometry, validate_rate,
)
from fasta import extract_id
from utils import DNAToBytes


MAGIC = b"GECC\x02\r\n\x1a"
METADATA_LIMIT = 16 * 1024 * 1024
FORMATS = ("v2", "legacy-genome", "legacy-chromosome")


@dataclass(frozen=True)
class Group:
    payload: bytes
    header: dict


@dataclass(frozen=True)
class ArchiveRecord:
    name: str
    title: str
    topology: str
    groups: tuple
    digest: str = ""

    @property
    def length(self):
        return sum(group.header["genome_length"] for group in self.groups)

    @property
    def blocks(self):
        return sum(group.header["block_num"] for group in self.groups)


@dataclass(frozen=True)
class Archive:
    format: str
    records: tuple


class RecoveryError(ValueError):
    """The supplied blocks could not be recovered and verified."""


@dataclass(frozen=True)
class RecoveredRecord:
    name: str
    title: str
    sequence: str
    missing_blocks: int
    rejected_blocks: int
    corrected_symbols: int


def _geometry(header):
    return tuple(header[key] for key in (
        "block_len", "chunk_size", "chunk_anchor_byte_num", "block_crc_byte_num",
    ))


def _json_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise EccFormatError(f"duplicate archive metadata field {key!r}")
        result[key] = value
    return result


def _validate_record_identity(name, title, topology, names):
    if (
        type(name) is not str or not name or name in names
        or type(title) is not str or any(char in title for char in "\r\n\x00")
        or extract_id(title) != name
        or topology not in ("linear", "circular", "unknown")
    ):
        raise EccFormatError(f"invalid or duplicate v2 identity for record {name!r}")


def _parse_v2(blob):
    start = len(MAGIC) + 4
    if not blob.startswith(MAGIC) or len(blob) < start + 32:
        raise EccFormatError("invalid or truncated GenomeECC v2 header")
    size = int.from_bytes(blob[len(MAGIC):start], "big")
    limit = len(blob) - 32
    if not 0 < size <= METADATA_LIMIT or start + size > limit:
        raise EccFormatError("invalid v2 metadata length")
    if sha256(memoryview(blob)[:limit]).digest() != blob[limit:]:
        raise EccFormatError("v2 archive checksum mismatch")
    try:
        metadata = json.loads(blob[start:start + size], object_pairs_hook=_json_object)
    except (UnicodeDecodeError, ValueError, RecursionError) as exc:
        raise EccFormatError(f"invalid v2 metadata: {exc}") from exc
    if type(metadata) is not list or not metadata:
        raise EccFormatError("v2 archive must contain at least one record")

    offset = start + size
    records, names = [], set()
    for item in metadata:
        if type(item) is not dict or item.keys() != {"id", "title", "topology", "sha256", "groups"}:
            raise EccFormatError("invalid v2 record metadata")
        name, title, topology, digest, sizes = (
            item["id"], item["title"], item["topology"], item["sha256"], item["groups"]
        )
        _validate_record_identity(name, title, topology, names)
        if (
            type(digest) is not str or len(digest) != 64
            or any(char not in "0123456789abcdef" for char in digest)
            or type(sizes) is not list or not sizes
        ):
            raise EccFormatError(f"invalid v2 metadata for record {name!r}")
        groups = []
        for i, length in enumerate(sizes):
            if type(length) is not int or not 0 < length <= limit - offset:
                raise EccFormatError(f"record {name!r}, group {i}: invalid payload length")
            payload = blob[offset:offset + length]
            offset += length
            header = unpack_chromosome(payload)
            if groups and _geometry(groups[0].header) != _geometry(header):
                raise EccFormatError(f"record {name!r}: inconsistent group geometry")
            if i < len(sizes) - 1 and header["genome_length"] % (
                header["block_len"] * header["chunk_size"]
            ):
                raise EccFormatError(f"record {name!r}: nonterminal group ends inside a chunk")
            groups.append(Group(payload, header))
        names.add(name)
        records.append(ArchiveRecord(name, title, topology, tuple(groups), digest))
    if offset != limit:
        raise EccFormatError("trailing bytes in v2 archive")
    return Archive("v2", tuple(records))


def parse_archive(blob, format="auto"):
    if format not in ("auto", *FORMATS):
        raise ValueError(f"unknown ECC format {format!r}")
    if format == "v2" or format == "auto" and blob.startswith(b"GECC"):
        return _parse_v2(blob)

    matches, errors = [], []
    for legacy in ("legacy-genome", "legacy-chromosome"):
        if format not in ("auto", legacy):
            continue
        try:
            if legacy == "legacy-genome":
                groups = unpack_legacy_genome(blob)
            else:
                groups = [(blob, unpack_chromosome(blob))]
        except EccFormatError as exc:
            errors.append(f"{legacy}: {exc}")
            continue
        records = tuple(
            ArchiveRecord(str(i), f"record-{i}", "unknown", (Group(payload, header),))
            for i, (payload, header) in enumerate(groups)
        )
        matches.append(Archive(legacy, records))
    if len(matches) == 1:
        return matches[0]
    if matches:
        raise EccFormatError("ambiguous legacy layout; select --format explicitly")
    raise EccFormatError("; ".join(errors))


def _sequence_digest(sequence):
    return sha256(sequence.upper().encode("ascii")).hexdigest()


def _encode_group(sequence, block_len, chunk_size, rate, mode, anchor, *, parity=None):
    chromosome = ChromosomeECC(block_len, chunk_size, rate, mode, anchor)
    chromosome.set_chr_seq(sequence)
    if parity is not None and parity != chromosome.redundant_block_num:
        chromosome.set_redundant_block_num(parity)
    return chromosome.get_all_ecc_bytes()


def encode_archive(records, *, block_len=320, chunk_size=10, rate=.04, mode="crc8",
                   anchor=8, format="v2", group_blocks=None):
    if format not in FORMATS:
        raise ValueError(f"unknown ECC format {format!r}")
    validate_geometry(block_len, chunk_size, anchor)
    validate_rate(rate)
    if group_blocks is not None and format != "v2":
        raise ValueError("group size is only supported by the v2 archive")
    if format == "v2":
        group_blocks = 4096 if group_blocks is None else group_blocks
        if type(group_blocks) is not int or not chunk_size <= group_blocks <= MAX_CODEWORD:
            raise ValueError(f"group size must be in {chunk_size}..{MAX_CODEWORD} blocks")
        group_blocks -= group_blocks % chunk_size
        validate_codeword(group_blocks, math.ceil(rate * group_blocks))
    metadata, payloads, names = [], [], set()
    for record in records:
        seq = record.sequence
        if not seq:
            raise ValueError(f"record {record.name!r}: empty sequence")
        if format == "legacy-chromosome" and payloads:
            raise ValueError("legacy chromosome format accepts exactly one record")
        if format == "legacy-genome" and len(payloads) == 255:
            raise ValueError("legacy genome format accepts at most 255 records")
        if format == "v2":
            _validate_record_identity(record.name, record.title, record.topology, names)
            names.add(record.name)
        step = group_blocks * block_len if format == "v2" else len(seq)
        sizes = []
        for start in range(0, len(seq), step):
            try:
                payload = _encode_group(seq[start:start + step], block_len, chunk_size, rate, mode, anchor)
            except ValueError as exc:
                raise ValueError(f"record {record.name!r}: {exc}") from exc
            sizes.append(len(payload))
            if format == "legacy-genome":
                payloads.append(len(payload).to_bytes(4, "big") + payload)
            else:
                payloads.append(payload)
        if format == "v2":
            metadata.append({
                "id": record.name, "title": record.title, "topology": record.topology,
                "sha256": _sequence_digest(seq), "groups": sizes,
            })
    if not payloads:
        raise ValueError("no records to encode")
    if format == "legacy-chromosome":
        return payloads[0]
    if format == "legacy-genome":
        return b''.join([bytes([len(payloads)]), *payloads])
    encoded_metadata = json.dumps(metadata, ensure_ascii=True, separators=(",", ":")).encode("ascii")
    if len(encoded_metadata) > METADATA_LIMIT:
        raise ValueError("archive metadata exceeds 16 MiB")
    blob = b''.join([MAGIC, len(encoded_metadata).to_bytes(4, "big"), encoded_metadata, *payloads])
    return blob + sha256(blob).digest()


def _matches_md5(sequence, header):
    packed = DNAToBytes(sequence + 'A' * (-len(sequence) % 4))
    expected = header["genome_md5"]
    if md5(packed).digest() == expected:
        return True
    # Pre-0.2.0 padded an already aligned sequence with an additional AAAA.
    return len(sequence) % 4 == 0 and md5(packed + b'\0').digest() == expected


def verify_records(archive, fasta_records, *, exact=False, block_len=None, chunk_size=None,
                   rate=None, mode=None, anchor=None):
    if rate is not None:
        validate_rate(rate)
    by_id = {record.name: record for record in archive.records} if archive.format == "v2" else None
    seen = set()
    for i, reference in enumerate(fasta_records):
        if i >= len(archive.records):
            raise ValueError("reference record count differs from archive")
        if by_id is None:
            record = archive.records[i]
        else:
            record = by_id.get(reference.name)
            if record is None or reference.name in seen:
                raise ValueError(f"unexpected or duplicate reference identifier {reference.name!r}")
        seen.add(record.name)
        if len(reference.sequence) != record.length:
            raise ValueError(f"record {record.name!r}: reference length mismatch")
        if record.digest and _sequence_digest(reference.sequence) != record.digest:
            raise ValueError(f"record {record.name!r}: SHA-256 mismatch")
        start = 0
        for group in record.groups:
            header = group.header
            block, chunk, anchor_bytes, crc = _geometry(header)
            crc_mode = CRC_MODES[crc]
            parity = header["redundant_block_num"]
            for key, actual, expected in (
                ("block length", block, block_len), ("chunk size", chunk, chunk_size),
                ("anchor width", anchor_bytes, anchor), ("CRC mode", crc_mode, mode),
            ):
                if expected is not None and actual != expected:
                    raise ValueError(f"record {record.name!r}: {key} differs from archive")
            if rate is not None and math.ceil(header["block_num"] * rate) != parity:
                raise ValueError(f"record {record.name!r}: redundancy differs from archive")
            end = start + header["genome_length"]
            sequence = reference.sequence[start:end]
            if not _matches_md5(sequence, header):
                raise ValueError(f"record {record.name!r}: sequence checksum mismatch")
            if exact and _encode_group(
                sequence, block, chunk, 0, crc_mode, anchor_bytes, parity=parity,
            ) != group.payload:
                raise ValueError(f"record {record.name!r}: regenerated ECC differs")
            start = end
    if len(seen) != len(archive.records):
        raise ValueError("reference record count differs from archive")


def _recover_group(group, candidates, offset):
    header = group.header
    chromosome = ChromosomeECC()
    chromosome.reset_with_unpacked_ecc(header)
    chromosome.rebuild_chunk_objs_with_ecc_codes(header["all_chunk_ecc_bytes"])
    block_len = chromosome.block_len
    length = chromosome.DNA_seq_len
    blank = 'A' * block_len
    sequences, erased = [], []
    missing = rejected = 0
    for i, expected_crc in enumerate(chromosome.block_crc_codes):
        index = offset + i
        sequence = candidates.get(index)
        if sequence is None:
            missing += 1
        else:
            valid = min(block_len, length - i * block_len)
            if len(sequence) not in (valid, block_len):
                raise RecoveryError(f"block {index}: expected {valid} bases")
            sequence = sequence[:valid].ljust(block_len, 'A')
            try:
                packed = DNAToBytes(sequence)
            except ValueError as exc:
                raise RecoveryError(f"block {index}: {exc}") from exc
            if block_crc(packed, chromosome.block_ecc_mode) == expected_crc:
                sequences.append(sequence)
                continue
            rejected += 1
        erased.append(i)
        sequences.append(blank)
    chromosome.block_seqs = sequences
    chromosome.block_seqs_to_bytes()
    try:
        corrected = chromosome.restore_colum_bytes_rs(header["all_rs_code_bytes"], erased)
    except ReedSolomonError as exc:
        raise RecoveryError(f"Reed–Solomon recovery failed: {exc}") from exc
    chromosome.colum_bytes_to_block_bytes()
    sequence = chromosome.assembled_chr_seq()
    if not _matches_md5(sequence, header):
        raise RecoveryError("recovered sequence checksum mismatch")
    return sequence, missing, rejected, corrected


def recover_blocks(archive, candidates):
    """Recover from a mapping of record IDs to zero-based block indices and DNA strings."""
    unknown = candidates.keys() - {record.name for record in archive.records}
    if unknown:
        raise RecoveryError(f"unknown record {next(iter(unknown))!r}")
    recovered = []
    for record in archive.records:
        blocks = candidates.get(record.name, {})
        if any(type(index) is not int or not 0 <= index < record.blocks for index in blocks):
            raise RecoveryError(f"record {record.name!r}: block index outside the record")
        parts = []
        offset = missing = rejected = corrected = 0
        for i, group in enumerate(record.groups):
            try:
                sequence, absent, bad, symbols = _recover_group(group, blocks, offset)
            except RecoveryError as exc:
                raise RecoveryError(f"record {record.name!r}, group {i}: {exc}") from exc
            parts.append(sequence)
            missing += absent
            rejected += bad
            corrected += symbols
            offset += group.header["block_num"]
        sequence = ''.join(parts)
        if record.digest and _sequence_digest(sequence) != record.digest:
            raise RecoveryError(f"record {record.name!r}: recovered SHA-256 mismatch")
        recovered.append(RecoveredRecord(record.name, record.title, sequence, missing, rejected, corrected))
    return recovered

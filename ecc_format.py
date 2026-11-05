"""Geometry and strict parsing for the historical chromosome and genome layouts."""

import math
import os
from pathlib import Path
import struct
from tempfile import NamedTemporaryFile


CHROMOSOME_HEADER = struct.Struct(">I16sHBH")
MAX_CODEWORD = 65_535
RS_POLYNOMIAL = 0x1002D
CRC_MODES = {1: "crc8", 2: "crc16"}
CRC_WIDTHS = {mode: width for width, mode in CRC_MODES.items()}


class EccFormatError(ValueError):
    """An ECC artifact has an invalid or unsupported layout."""


def ceil_div(value, divisor):
    return (value + divisor - 1) // divisor


def validate_geometry(block_len, chunk_size, anchor_bytes):
    if type(block_len) is not int or not 4 <= block_len <= 1020 or block_len % 4:
        raise ValueError("block length must be a multiple of 4 in 4..1020")
    if type(chunk_size) is not int or not 1 <= chunk_size <= 63:
        raise ValueError("chunk size must be in 1..63")
    if type(anchor_bytes) is not int or not 0 <= anchor_bytes <= 15:
        raise ValueError("anchor width must be in 0..15 bytes")
    if anchor_bytes > (block_len // 4) * chunk_size:
        raise ValueError("anchor width exceeds the chunk capacity")


def validate_rate(rate):
    if not math.isfinite(rate) or not 0 <= rate <= 1:
        raise ValueError("redundancy must be finite and in 0..1")


def validate_codeword(blocks, parity):
    if type(parity) is not int or not 0 <= parity < MAX_CODEWORD:
        raise ValueError("invalid Reed–Solomon parity symbol count")
    if blocks + parity > MAX_CODEWORD:
        raise ValueError(
            f"Reed–Solomon codeword has {blocks} data + {parity} parity symbols; "
            f"maximum is {MAX_CODEWORD}. Use a grouped archive."
        )


def unpack_chromosome(blob):
    if len(blob) < CHROMOSOME_HEADER.size:
        raise EccFormatError("truncated chromosome ECC header")
    length, digest, geometry, widths, parity = CHROMOSOME_HEADER.unpack_from(blob)
    block_len, chunk_size = geometry >> 6, geometry & 63
    anchor_bytes, crc_bytes = widths >> 4, widths & 15
    try:
        validate_geometry(block_len, chunk_size, anchor_bytes)
        if not length:
            raise ValueError("chromosome length must be positive")
        if crc_bytes not in CRC_MODES:
            raise ValueError(f"unsupported block CRC width {crc_bytes}")
        blocks = ceil_div(length, block_len)
        validate_codeword(blocks, parity)
    except ValueError as exc:
        raise EccFormatError(str(exc)) from exc
    chunks = ceil_div(blocks, chunk_size)
    block_size = block_len // 4
    chunk_end = CHROMOSOME_HEADER.size + chunks * (anchor_bytes + chunk_size * crc_bytes)
    expected = chunk_end + 2 * ceil_div(block_size, 2) * parity
    if len(blob) != expected:
        raise EccFormatError(
            f"chromosome ECC length mismatch: expected {expected}, got {len(blob)}"
        )
    return {
        "genome_length": length,
        "genome_md5": bytes(digest),
        "block_len": block_len,
        "chunk_size": chunk_size,
        "chunk_anchor_byte_num": anchor_bytes,
        "block_crc_byte_num": crc_bytes,
        "block_size": block_size,
        "redundant_block_num": parity,
        "block_num": blocks,
        "chunk_num": chunks,
        "all_chunk_ecc_bytes": bytes(blob[CHROMOSOME_HEADER.size:chunk_end]),
        "all_rs_code_bytes": bytes(blob[chunk_end:]),
        "chunk_objs": [],
    }


def unpack_legacy_genome(blob):
    if not blob or not blob[0]:
        raise EccFormatError("legacy genome ECC must contain at least one record")
    records = []
    offset = 1
    for index in range(blob[0]):
        if len(blob) - offset < 4:
            raise EccFormatError(f"record {index}: truncated ECC length field")
        size = int.from_bytes(blob[offset:offset + 4], "big")
        offset += 4
        if size > len(blob) - offset:
            raise EccFormatError(f"record {index}: truncated chromosome ECC payload")
        payload = blob[offset:offset + size]
        try:
            header = unpack_chromosome(payload)
        except EccFormatError as exc:
            raise EccFormatError(f"record {index}: {exc}") from exc
        records.append((payload, header))
        offset += size
    if offset != len(blob):
        raise EccFormatError("trailing bytes in legacy genome ECC")
    return records


def write_atomic(path, data):
    """Replace the destination only after writing the complete binary result."""
    path = Path(path)
    temporary = None
    try:
        with NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", delete=False) as file:
            temporary = Path(file.name)
            if isinstance(data, bytes):
                file.write(data)
            else:
                file.writelines(data)
            file.flush()
            os.fsync(file.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)

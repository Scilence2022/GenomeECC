"""Archive compatibility and recovery through public API and CLI."""

from hashlib import md5
from pathlib import Path
from random import Random
import subprocess
import sys

import pytest

from chromosomeChunk import block_crc
from chromosomeECC import ChromosomeECC
from ecc_archive import RecoveryError, encode_archive, parse_archive, recover_blocks, verify_records
from ecc_format import EccFormatError
from fasta import FastaRecord, read_fasta
from genomeECC import GenomeECC
from utils import DNAToBytes, bytesToDNA


ROOT = Path(__file__).resolve().parents[1]


def record(name, length, seed=7, topology="linear"):
    sequence = ''.join(Random(seed).choices("ACGT", k=length))
    return FastaRecord(name, f"{name} topology={topology}", sequence, topology)


def candidate_blocks(sequence, width):
    return {i // width: sequence[i:i + width] for i in range(0, len(sequence), width)}


def cli(*args):
    return subprocess.run(
        [sys.executable, str(ROOT / "genome_ecc_cli.py"), *map(str, args)],
        capture_output=True, text=True,
    )


@pytest.mark.parametrize("format", ["v2", "legacy-genome", "legacy-chromosome"])
def test_roundtrip_formats_and_checksum(format):
    source = record("alpha", 321)
    data = encode_archive([source], format=format, mode="crc16", anchor=6)
    archive = parse_archive(data)
    assert archive.format == format
    assert archive.records[0].length == 321
    verify_records(archive, [source], exact=True)
    # The entire last block contains one real base.
    blocks = candidate_blocks(source.sequence, 320)
    del blocks[1]
    recovered, = recover_blocks(archive, {archive.records[0].name: blocks})
    assert recovered.sequence == source.sequence
    assert recovered.missing_blocks == 1


def test_v2_preserves_identity_topology_and_group_geometry():
    alpha, beta = record("alpha", 101, topology="circular"), record("beta", 40, seed=6)
    archive = parse_archive(encode_archive(
        [alpha, beta], block_len=12, chunk_size=2, anchor=4, group_blocks=5, rate=.5,
    ))
    assert [r.name for r in archive.records] == ["alpha", "beta"]
    assert archive.records[0].topology == "circular"
    assert [g.header["block_num"] for g in archive.records[0].groups] == [4, 4, 1]
    verify_records(archive, [beta, alpha], exact=True)
    blocks = candidate_blocks(alpha.sequence, 12)
    for index in (0, 3, 5, 8):
        del blocks[index]
    restored = recover_blocks(archive, {"alpha": blocks, "beta": candidate_blocks(beta.sequence, 12)})
    assert [r.sequence for r in restored] == [alpha.sequence, beta.sequence]
    assert restored[0].missing_blocks == 4

    del blocks[1]
    with pytest.raises(RecoveryError, match="group 0.*too many erasures"):
        recover_blocks(archive, {"alpha": blocks})


def test_codeword_limit_and_grouping():
    source = record("long", 65_535 * 4)
    config = dict(block_len=4, chunk_size=63, anchor=1, rate=1 / 65_535)
    with pytest.raises(ValueError, match="maximum is 65535"):
        encode_archive([source], format="legacy-chromosome", **config)
    grouped = parse_archive(encode_archive([source], **config))
    assert grouped.records[0].length == len(source.sequence)
    assert len(grouped.records[0].groups) > 1
    verify_records(grouped, [source])


def test_last_legal_single_codeword_corrects_an_erasure():
    source = record("limit", 65_534 * 4)
    chromosome = ChromosomeECC(4, 63, 1 / 65_535, chunk_anchor_byte_num=1)
    chromosome.set_chr_seq(source.sequence)
    parity = chromosome.get_colum_rs(0)
    column = chromosome.colum(0)
    chromosome.colum_bytes[0] = b'\0\0' + column[2:]
    chromosome.restore_colum_bytes_rs(parity, [0])
    assert chromosome.colum(0) == column
    assert len(parity) == 2


def test_short_anchor_is_padded_and_zero_redundancy_passes_clean_blocks():
    source = record("x", 4)
    archive = parse_archive(encode_archive([source], rate=0))
    group, = archive.records[0].groups
    header = group.header
    assert header["all_chunk_ecc_bytes"][:8] == DNAToBytes(source.sequence) + b'\0' * 7
    assert header["all_rs_code_bytes"] == b''
    result, = recover_blocks(archive, {"x": {0: source.sequence}})
    assert result.sequence == source.sequence
    with pytest.raises(RecoveryError, match="too many erasures"):
        recover_blocks(archive, {})


def test_recovery_corrects_mixed_erasure_and_crc_collision():
    source = record("x", 48)
    archive = parse_archive(encode_archive(
        [source], block_len=8, chunk_size=2, anchor=1, rate=.5,
    ))
    blocks = candidate_blocks(source.sequence, 8)
    expected = DNAToBytes(blocks[0])
    crc = block_crc(expected, "crc8")
    collision = next(
        value for n in range(65_536) if (value := n.to_bytes(2, "big")) != expected
        and block_crc(value, "crc8") == crc
    )
    blocks[0] = bytesToDNA(collision)
    del blocks[5]
    result, = recover_blocks(archive, {"x": blocks})
    assert result.sequence == source.sequence
    assert result.missing_blocks == 1
    assert result.rejected_blocks == 0
    assert result.corrected_symbols == 2
    del blocks[4]
    with pytest.raises(RecoveryError):
        recover_blocks(archive, {"x": blocks})


def test_crc_rejections_count_as_erasures():
    source = record("x", 800)
    archive = parse_archive(encode_archive([source], mode="crc16"))
    blocks = candidate_blocks(source.sequence, 320)
    blocks[0] = ('T' if blocks[0][0] != 'T' else 'A') + blocks[0][1:]
    restored, = recover_blocks(archive, {"x": blocks})
    assert restored.sequence == source.sequence
    assert restored.rejected_blocks == 1


@pytest.mark.parametrize("format", ["v2", "legacy-genome", "legacy-chromosome"])
def test_reader_rejects_truncated_and_extra_bytes(format):
    blob = encode_archive([record("alpha", 700)], format=format)
    for data in (b'', b'\0ignored', blob[:4], blob[:-1], blob + b'ignored'):
        with pytest.raises(EccFormatError):
            parse_archive(data, format)
    if format == "v2":
        corrupted = bytearray(blob)
        corrupted[-33] ^= 1
        with pytest.raises(EccFormatError, match="checksum mismatch"):
            parse_archive(bytes(corrupted))


def test_reader_detects_repository_legacy_samples():
    samples = sorted(p for p in (ROOT / "genomes").glob("*") if p.suffix.lower() == ".ecc")
    if not samples:
        pytest.skip("repository fixtures are absent from source distributions")
    assert len(samples) == 5
    for sample in samples:
        archive = parse_archive(sample.read_bytes())
        assert archive.format == "legacy-chromosome"
        assert len(archive.records) == 1
    paired = parse_archive((ROOT / "genomes/cgl-ssbl-from_nxm.fna.ecc").read_bytes())
    verify_records(paired, read_fasta(ROOT / "genomes/cgl-ssbl-from_nxm.fna"))


def test_legacy_aligned_md5_can_be_verified_and_recovered():
    source = record("x", 640)
    data = encode_archive([source], format="legacy-chromosome")
    historical = data[:4] + md5(DNAToBytes(source.sequence) + b'\0').digest() + data[20:]
    archive = parse_archive(historical)
    verify_records(archive, [source])
    with pytest.raises(ValueError, match="regenerated ECC differs"):
        verify_records(archive, [source], exact=True)
    result, = recover_blocks(archive, {"0": candidate_blocks(source.sequence, 320)})
    assert result.sequence == source.sequence


def test_configuration_change_rebuilds_derived_data():
    source = record("x", 321)
    chromosome = ChromosomeECC(320, 10)
    chromosome.set_chr_seq(source.sequence)
    chromosome.get_all_ecc_bytes()
    chromosome.set_block_ecc_mode("crc16")
    chromosome.set_redundant_block_num(2)
    chromosome.set_chunk_size(1)
    archive = parse_archive(chromosome.get_all_ecc_bytes())
    header = archive.records[0].groups[0].header
    assert (header["block_crc_byte_num"], header["chunk_num"], header["redundant_block_num"]) == (2, 2, 2)
    verify_records(archive, [source], exact=True)


@pytest.mark.parametrize("body, error", [
    ("AAAA\n>x\nAAAA\n", "before first FASTA header"),
    (">x\n>x2\nAAAA\n", "empty FASTA record"),
    (">x\nAAAA\n>x\nAAAA\n", "duplicate FASTA identifier"),
    (">x\nACGT\nAN\n", "base 6"),
])
def test_fasta_reports_malformed_input(tmp_path, body, error):
    path = tmp_path / "bad.fa"
    path.write_text(body)
    with pytest.raises(ValueError, match=error):
        list(read_fasta(path))


def test_reusing_genome_replaces_old_records_and_cache(tmp_path):
    path = tmp_path / "genome.fa"
    genome = GenomeECC()
    genome.add_chromosome("old", "ACGT" * 200)
    genome.genChromesomeECCs()
    path.write_text(">new non-circular\naaaa\n")
    genome.read_genome(path)
    genome.genChromesomeECCs()
    assert list(genome.chr_ECC_objs) == ["new"]
    assert not genome.chr_ECC_objs["new"].circular_DNA
    assert len(parse_archive(genome.get_all_genome_ecc()).records) == 1
    with pytest.raises(ValueError, match="duplicate"):
        genome.read_chromosome(path)


def test_cli_output_is_atomic_and_recovers_without_reference(tmp_path):
    source = record("x", 37)
    reference, ecc, blocks, output = [tmp_path / name for name in ("ref.fa", "test.ecc", "blocks.tsv", "out.fa")]
    reference.write_text(f">{source.title}\n{source.sequence}\n")
    assert cli("encode", "-g", reference, "-o", ecc, "--block-len", 12, "--chunk-size", 2,
               "--anchor-bytes", 2, "--redundancy", ".25").returncode == 0
    assert cli("verify", "-g", reference, "-e", ecc, "--exact").returncode == 0
    assert cli("info", "-e", ecc).returncode == 0
    reference.unlink()
    candidates = candidate_blocks(source.sequence, 12)
    output.write_text("existing output")
    blocks.write_text("")
    failure = cli("recover", "-e", ecc, "--blocks", blocks, "-o", output)
    assert failure.returncode == 1
    assert "too many erasures" in failure.stderr
    assert output.read_text() == "existing output"

    blocks.write_text(''.join(f"x\t{i}\t{value}\n" for i, value in candidates.items() if i != 3))
    recovered = cli("recover", "-e", ecc, "--blocks", blocks, "-o", output)
    assert recovered.returncode == 0, recovered.stderr
    assert list(read_fasta(output)) == [source]


def test_cli_migrates_only_with_a_matching_reference(tmp_path):
    source = record("x", 37)
    ecc, reference, output = [tmp_path / name for name in ("old.ecc", "ref.fa", "v2.ecc")]
    ecc.write_bytes(encode_archive([source], format="legacy-chromosome"))
    reference.write_text(f">{source.title}\n{'A' * 37}\n")
    rejected = cli("migrate", "-e", ecc, "-g", reference, "-o", output)
    assert rejected.returncode == 1
    assert not output.exists()

    reference.write_text(f">{source.title}\n{source.sequence}\n")
    accepted = cli("migrate", "-e", ecc, "-g", reference, "-o", output)
    assert accepted.returncode == 0, accepted.stderr
    assert parse_archive(output.read_bytes()).records[0].name == "x"

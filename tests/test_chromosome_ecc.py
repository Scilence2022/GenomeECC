"""Core ECC invariants for a single chromosome.

The central property: a genome must survive block decomposition, column
decomposition, Reed-Solomon encoding, damage, and reassembly unchanged.
"""
import random

import pytest

from chromosomeECC import ChromosomeECC


def random_dna(n, seed=0):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


def build(seq, block_len=320, chunk_size=6, rate=0.20):
    c = ChromosomeECC(block_len=block_len, chunk_size=chunk_size,
                      redundant_block_rate=rate)
    c.set_chr_seq(seq)
    c.block_bytes_to_colum_bytes()
    return c


def columns_from_blocks(c):
    """Rebuild columns from (possibly damaged) block_bytes."""
    n = c.colum_byte_num_per_block
    c.colum_bytes = [
        b''.join(c._padded_block_bytes(j)[i * n:(i + 1) * n] for j in range(c.block_num))
        for i in range(c.colum_num())
    ]


# block_len=500 and 300 give an ODD block_size (125, 75). Those are the defaults of
# ChromosomeECC and GenomeECC respectively, and used to lose the last byte of every
# block: it belonged to no column, so it was neither RS-protected nor restored.
BLOCK_LENS = [160, 200, 300, 320, 500, 520]


@pytest.mark.parametrize("block_len", BLOCK_LENS)
def test_columns_cover_every_byte_of_every_block(block_len):
    c = build(random_dna(20000, 1), block_len=block_len)
    covered = c.colum_num() * c.colum_byte_num_per_block
    assert covered >= c.block_size, (
        "columns cover %d of %d bytes per block -- the remainder is unprotected"
        % (covered, c.block_size))


@pytest.mark.parametrize("block_len", BLOCK_LENS)
def test_block_column_round_trip_is_lossless(block_len):
    seq = random_dna(20000, 2)
    c = build(seq, block_len=block_len)
    original = list(c.block_bytes)
    c.colum_bytes_to_block_bytes()
    assert c.block_bytes == original


@pytest.mark.parametrize("block_len", BLOCK_LENS)
def test_assembled_sequence_matches_original(block_len):
    seq = random_dna(20000, 3)
    c = build(seq, block_len=block_len)
    assert c.assembled_chr_seq() == seq


@pytest.mark.parametrize("block_len", BLOCK_LENS)
def test_rs_recovers_wiped_blocks(block_len):
    """The core value proposition: erased blocks are reconstructed exactly."""
    seq = random_dna(20000, 4)
    ref = build(seq, block_len=block_len)
    rs = b''.join(ref.get_colum_rs(i) for i in range(ref.colum_num()))

    damaged = build(seq, block_len=block_len)
    dead = sorted(random.Random(block_len).sample(range(damaged.block_num), 5))
    for b in dead:
        damaged.block_bytes[b] = bytes(damaged.block_size)
    columns_from_blocks(damaged)

    damaged.restore_colum_bytes_rs(rs, erasure_positions=dead)
    damaged.colum_bytes_to_block_bytes()
    assert damaged.assembled_chr_seq() == seq


def test_rs_corrects_up_to_nsym_erasures_and_no_more():
    seq = random_dna(20000, 5)
    ref = build(seq, block_len=520, rate=0.20)
    nsym = ref.redundant_block_num
    rs = b''.join(ref.get_colum_rs(i) for i in range(ref.colum_num()))

    def attempt(n_dead):
        d = build(seq, block_len=520, rate=0.20)
        dead = sorted(random.Random(n_dead).sample(range(d.block_num), n_dead))
        for b in dead:
            d.block_bytes[b] = bytes(d.block_size)
        columns_from_blocks(d)
        d.restore_colum_bytes_rs(rs, erasure_positions=dead)
        d.colum_bytes_to_block_bytes()
        return d.assembled_chr_seq() == seq

    assert attempt(nsym) is True
    with pytest.raises(Exception):
        attempt(nsym + 1)


class TestHeaderPacking:
    @pytest.mark.parametrize("block_len,chunk_size", [(500, 6), (300, 10), (1023, 63), (4, 1)])
    def test_pack_unpack_round_trip(self, block_len, chunk_size):
        c = ChromosomeECC(block_len=block_len, chunk_size=chunk_size)
        packed = c.pack_block_len_chunk_size()
        assert c.unpack_block_len_chunk_size(packed) == [block_len, chunk_size]

    def test_rejects_block_len_over_10_bits(self):
        c = ChromosomeECC(block_len=1024, chunk_size=6)
        with pytest.raises(ValueError, match="10 bits"):
            c.pack_block_len_chunk_size()

    def test_rejects_chunk_size_over_6_bits(self):
        # 500/64 silently unpacked as 501/0 before validation existed.
        c = ChromosomeECC(block_len=500, chunk_size=64)
        with pytest.raises(ValueError, match="6 bits"):
            c.pack_block_len_chunk_size()


class TestEccModes:
    def test_crc8_declares_one_byte(self):
        # The default mode is crc8 but block_ecc_byte_num was hardcoded to 2, so the
        # header advertised a 2-byte CRC width for 1-byte digests.
        c = ChromosomeECC()
        assert c.block_ecc_mode == 'crc8'
        assert c.block_ecc_byte_num == 1

    @pytest.mark.parametrize("mode,width", [('crc8', 1), ('crc16', 2)])
    def test_digest_width_matches_declared_width(self, mode, width):
        c = ChromosomeECC(block_len=320, chunk_size=6, block_ecc_mode=mode)
        c.set_chr_seq(random_dna(4000, 6))
        assert c.block_ecc_byte_num == width
        assert all(len(code) == width for code in c.block_crc_codes)

    def test_unknown_mode_raises(self):
        c = ChromosomeECC(block_len=320)
        c.block_ecc_mode = 'crc99'
        with pytest.raises(ValueError, match="unknown block_ecc_mode"):
            c.set_chr_seq(random_dna(1000, 7))


def test_round_dna_pads_only_when_needed():
    c = ChromosomeECC()
    assert len(c.round_DNA('ACGT')) == 4        # already aligned: no padding
    assert len(c.round_DNA('ACG')) == 4
    assert len(c.round_DNA('A')) == 4
    assert c.round_DNA('ACGT') == 'ACGT'


def test_set_chr_seq_is_idempotent():
    c = ChromosomeECC(block_len=320, chunk_size=6)
    seq = random_dna(8000, 8)
    c.set_chr_seq(seq)
    crcs, anchors = list(c.block_crc_codes), list(c.chunk_anchores)
    c.set_chr_seq(seq)
    assert c.block_crc_codes == crcs
    assert c.chunk_anchores == anchors
    assert len(c.block_crc_codes) == c.block_num


def test_anchor_width_follows_chunk_anchor_byte_num():
    c = ChromosomeECC(block_len=320, chunk_size=6, chunk_anchor_byte_num=8)
    c.set_chr_seq(random_dna(8000, 9))
    assert all(len(a) == 8 for a in c.chunk_anchores)
    assert len(c.get_chunk_anchor(0)) == 8


def test_chunk_block_crc_codes_pads_short_final_chunk():
    c = ChromosomeECC(block_len=320, chunk_size=6)
    c.set_chr_seq(random_dna(8000, 10))
    last = c.get_chunk_block_crc_codes(c.chunk_num - 1)
    assert len(last) == c.chunk_size
    assert all(len(code) == c.block_ecc_byte_num for code in last)


def test_zero_crc_does_not_trigger_recompute_loop():
    # b'\x00' is falsy; the accessors used to treat a valid zero CRC as "not computed".
    c = ChromosomeECC(block_len=320, chunk_size=6)
    c.set_chr_seq(random_dna(8000, 11))
    c.block_crc_codes[0] = bytes(c.block_ecc_byte_num)
    assert c.get_block_crc(0) == bytes(c.block_ecc_byte_num)
    assert len(c.block_crc_codes) == c.block_num


def test_rscodec_keeps_gf16_when_redundancy_changes():
    c = ChromosomeECC(block_len=320, chunk_size=6)
    c.set_chr_seq(random_dna(8000, 12))
    c.set_redundant_block_rate(0.05)
    assert c.ecc.nsize == c.nsize        # rebuilding used to drop it back to 255


def test_ecc_bytes_are_deterministic():
    seq = random_dna(8000, 13)
    a = ChromosomeECC(block_len=320, chunk_size=6)
    a.set_chr_seq(seq)
    b = ChromosomeECC(block_len=320, chunk_size=6)
    b.set_chr_seq(seq)
    assert a.get_all_ecc_bytes() == b.get_all_ecc_bytes()


def test_unpack_ecc_header_matches_encoder():
    c = ChromosomeECC(block_len=320, chunk_size=6, chunk_anchor_byte_num=8)
    seq = random_dna(8000, 14)
    c.set_chr_seq(seq)
    unpacked = c.unpack_ecc_bytes(c.get_all_ecc_bytes())
    assert unpacked['genome_length'] == len(seq)
    assert unpacked['block_len'] == 320
    assert unpacked['chunk_size'] == 6
    assert unpacked['chunk_anchor_byte_num'] == 8
    assert unpacked['block_crc_byte_num'] == c.block_ecc_byte_num
    assert unpacked['block_num'] == c.block_num
    assert unpacked['chunk_num'] == c.chunk_num


def test_md5_covers_only_the_real_sequence():
    # round_DNA used to append a spurious 4-base group for aligned lengths, which fed a
    # phantom byte into the stored MD5.
    from hashlib import md5
    from utils import DNAToBytes
    seq = random_dna(8000, 15)          # 8000 % 4 == 0
    c = ChromosomeECC(block_len=320, chunk_size=6)
    c.set_chr_seq(seq)
    assert c.get_chr_md5() == md5(DNAToBytes(seq)).digest()

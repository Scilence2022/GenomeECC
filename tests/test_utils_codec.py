"""DNA <-> bytes codec, XOR and byte-permutation invariants."""
import random

import pytest

from utils import (DNAToByte, DNAToBytes, byteToDNA, bytesToDNA, calc_gc,
                   check_dna, max_homo_len, rand_order,
                   reverse_sort_bytes_by_order, rev_seq, sha256_rand,
                   sort_bytes_by_order, xor)


def random_dna(n, seed=0):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


@pytest.mark.parametrize("n", [4, 8, 400, 4000])
def test_dna_bytes_round_trip(n):
    dna = random_dna(n, seed=n)
    assert bytesToDNA(DNAToBytes(dna)) == dna


def test_codec_bit_layout():
    # A=00 C=01 G=10 T=11, first base in the high bits.
    assert DNAToByte('ACGT') == bytes([0b00011011])
    assert byteToDNA(0b00011011) == 'ACGT'
    assert DNAToBytes('AAAA') == b'\x00'
    assert DNAToBytes('TTTT') == b'\xff'


def test_every_byte_round_trips():
    all_bytes = bytes(range(256))
    assert DNAToBytes(bytesToDNA(all_bytes)) == all_bytes


def test_dna_to_bytes_is_case_insensitive():
    assert DNAToBytes('acgt') == DNAToBytes('ACGT')


def test_dna_to_bytes_rejects_non_acgt():
    # Real genomes contain N; silently mapping it to a base would corrupt data.
    with pytest.raises(ValueError, match="non-ACGT"):
        DNAToBytes('ACGN')


def test_dna_to_bytes_rejects_unaligned_length():
    with pytest.raises(ValueError, match="multiple of 4"):
        DNAToBytes('ACG')


def test_dna_to_bytes_empty():
    assert DNAToBytes('') == b''
    assert bytesToDNA(b'') == ''


class TestXor:
    def test_equal_length(self):
        assert xor(bytes([0xff] * 8), bytes([0x0f] * 8)) == bytes([0xf0] * 8)

    def test_shorter_operand_is_right_zero_padded(self):
        assert xor(b'\x01\x02', b'\x03') == b'\x02\x02'
        assert xor(b'\x03', b'\x01\x02') == b'\x02\x02'

    def test_empty(self):
        assert xor(b'', b'') == b''

    def test_involution(self):
        rng = random.Random(3)
        a = bytes(rng.randrange(256) for _ in range(64))
        b = bytes(rng.randrange(256) for _ in range(64))
        assert xor(xor(a, b), b) == a


def test_sort_bytes_by_order_round_trip():
    data = bytes(range(32))
    order = rand_order(32, 7)
    shuffled = sort_bytes_by_order(data, order)
    assert shuffled != data
    assert reverse_sort_bytes_by_order(shuffled, order) == data


def test_sha256_rand_accepts_large_seeds():
    # Previously bytes([seed]) raised ValueError for anything above 255.
    assert sha256_rand(300) != sha256_rand(301)
    assert isinstance(sha256_rand(0), int)


def test_gc_and_homopolymer():
    assert calc_gc('GCGC') == 1.0
    assert calc_gc('ATAT') == 0.0
    assert calc_gc('ACGT') == 0.5
    assert max_homo_len('AAACGT') == 3
    assert max_homo_len('ACGTTTTT') == 5
    assert max_homo_len('ACGT') == 1


def test_check_dna_bounds():
    assert check_dna('ACGT' * 10)
    assert not check_dna('GC' * 20)          # GC content 1.0
    assert not check_dna('A' * 40)           # homopolymer + GC 0.0


def test_rev_seq_is_reverse_complement():
    assert rev_seq('AAAACCCGGT') == 'ACCGGGTTTT'
    assert rev_seq(rev_seq('ACGTACGTAC')) == 'ACGTACGTAC'

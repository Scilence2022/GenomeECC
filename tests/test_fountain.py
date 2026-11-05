"""Fountain-code degree distributions and the Glass decoder."""
import random

import numpy as np
import pytest

from fountain import Fountain
from glass import Glass
from utils import ideal_distribution, robust_distribution


class FakeDroplet:
    def __init__(self, chunk_nums, data):
        self._nums = list(chunk_nums)
        self.data = data

    def get_chunk_nums(self):
        return self._nums


class TestDistributions:
    @pytest.mark.parametrize("K", [2, 10, 50, 100, 500, 1000])
    def test_robust_distribution_is_valid_for_any_K(self, K):
        # Used to raise a broadcast ValueError whenever M > K, which the code's own
        # default c_value=0.01 triggers for every K below ~500.
        p = robust_distribution(K, 0.01, 0.01)
        assert len(p) == K
        assert np.isclose(np.sum(p), 1.0)
        assert np.all(p >= 0)

    @pytest.mark.parametrize("K", [10, 1000])
    @pytest.mark.parametrize("c,delta", [(0.01, 0.01), (0.2, 0.1), (0.5, 0.1)])
    def test_robust_distribution_across_parameters(self, K, c, delta):
        p = robust_distribution(K, c, delta)
        assert len(p) == K
        assert np.isclose(np.sum(p), 1.0)

    @pytest.mark.parametrize("K", [2, 10, 1000])
    def test_ideal_distribution_normalised(self, K):
        p = ideal_distribution(K)
        assert len(p) == K
        assert np.isclose(np.sum(p), 1.0)

    def test_distributions_reject_degenerate_K(self):
        with pytest.raises(ValueError):
            ideal_distribution(1)
        with pytest.raises(ValueError):
            robust_distribution(1)


class TestFountain:
    def test_constructs_for_small_chunk_counts(self):
        f = Fountain(data=b'x' * 350, chunk_size=35)   # K=10
        assert f.num_of_chunks == 10
        assert len(f.degrees) == 30

    def test_gen_degrees_arity(self):
        # gen_degrees passed 7 args to a 6-parameter function: unconditional TypeError.
        f = Fountain(data=b'x' * 35000, chunk_size=35)
        f.gen_degrees()
        assert len(f.degrees) == 3000


class TestGlassDecode:
    def test_peels_a_chain(self):
        g = Glass(3)
        g.addDroplet(FakeDroplet([0], bytes([0x01])))
        g.addDroplet(FakeDroplet([0, 1], bytes([0x01 ^ 0x02])))
        g.addDroplet(FakeDroplet([1, 2], bytes([0x02 ^ 0x04])))
        assert g.decode() is True
        assert g.chunks == [bytes([1]), bytes([2]), bytes([4])]

    def test_adjacent_singletons_are_not_skipped(self):
        # decode() popped from self.entries while enumerating it, skipping entries.
        g = Glass(3)
        for i, val in enumerate([b'\x01', b'\x02', b'\x03']):
            g.addDroplet(FakeDroplet([i], val))
        assert g.decode() is True
        assert g.chunks == [b'\x01', b'\x02', b'\x03']

    def test_decodes_with_few_droplets(self):
        # The old code refused outright at <= 10 droplets and left every chunk None.
        g = Glass(1)
        g.addDroplet(FakeDroplet([0], b'\xaa'))
        assert g.decode() is True
        assert g.chunks == [b'\xaa']

    def test_terminates_without_any_degree_one_droplet(self):
        # This spun forever before: no degree-1 entry means no progress, and the loop
        # condition only exited once something had been solved.
        g = Glass(2)
        for _ in range(12):
            g.addDroplet(FakeDroplet([0, 1], b'\x0f'))
        assert g.decode() is False
        assert g.chunks == [None, None]

    def test_does_not_mutate_droplet_state(self):
        d = FakeDroplet([0, 1], b'\x0f')
        g = Glass(2)
        g.addDroplet(d)
        g.addDroplet(FakeDroplet([0], b'\x0a'))
        g.decode()
        assert d._nums == [0, 1]

    def test_get_string_returns_bytes(self):
        g = Glass(2)
        g.addDroplet(FakeDroplet([0], b'\xaa'))
        g.decode()
        assert isinstance(g.getString(), bytes)

    def test_write_to_file_refuses_incomplete_decode(self, tmp_path):
        g = Glass(2)
        g.addDroplet(FakeDroplet([0], b'\xaa'))
        g.decode()
        assert not g.isDone()
        with pytest.raises(ValueError, match="unrecovered"):
            g.writeToFile(str(tmp_path / "out.bin"))

    def test_write_to_file_when_complete(self, tmp_path):
        g = Glass(2)
        g.addDroplet(FakeDroplet([0], b'\xaa'))
        g.addDroplet(FakeDroplet([1], b'\xbb'))
        g.decode()
        out = tmp_path / "out.bin"
        assert g.writeToFile(str(out)) is True
        assert out.read_bytes() == b'\xaa\xbb'


def test_importing_glass_does_not_change_recursion_limit():
    import sys
    assert sys.getrecursionlimit() < 100000


def test_incremental_peeling_and_conflict_detection():
    decoder = Glass(2)
    decoder.addDroplet(FakeDroplet([0], b'\x01'))
    assert not decoder.decode()
    decoder.addDroplet(FakeDroplet([0, 1], b'\x03'))
    assert decoder.decode()
    assert decoder.getString() == b'\x01\x02'
    with pytest.raises(ValueError, match="conflicting"):
        decoder.addDroplet(FakeDroplet([0], b'\x10'))
    assert not decoder.isDone()


def test_wire_checks_actual_index_without_mutating_template():
    from DNAdroplet import DNADroplet

    droplet = DNADroplet(b'abcd', 42)
    wire = droplet.to_DNA()
    template = DNADroplet(bytes(4), 42)
    assert not template.set_droplet_from_DNA('T' + wire[1:])
    assert not template.set_droplet_from_DNA(wire[:-4])
    assert template.data == bytes(4)
    template = DNADroplet(bytes(4))
    assert template.set_droplet_from_DNA(wire)
    assert template.index == 42
    assert template.data == b'abcd'


def test_fountain_legacy_schedule_and_wire_vector():
    from DNAfountain import DNAFountain

    fountain = DNAFountain(b'abcdefghijklmnopqrst', chunk_size=4, seed=11)
    assert fountain.degrees == [1, 1, 3, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2]
    assert fountain.DNAdroplet().to_DNA() == "AAAAAAAAAAAAAAAGCTTCCTTGCTTTCGAAATTCTATCGTCCCTAT"


@pytest.mark.parametrize("data", [b'', b'\0', b'ABCDEFG\0xyz\0\0'])
def test_manifest_fountain_round_trip_with_dropped_strands(data):
    from fountain_codec import FountainManifest, decode_dna, encode_dna

    manifest, strands = encode_dna(data, 100 if data else 0, chunk_size=3, seed=4)
    manifest = FountainManifest.from_bytes(manifest.to_bytes())
    assert decode_dna(manifest, strands[::2]) == data
    if data:
        with pytest.raises(ValueError, match="unrecovered"):
            decode_dna(manifest, [])


def test_fountain_filter_has_attempt_limit():
    from fountain_codec import encode_dna

    with pytest.raises(ValueError, match="budget exhausted"):
        encode_dna(b'abc', 1, strand_filter=lambda _: False, max_attempts=2)


def test_soliton_bounds_and_random_state():
    from utils import get_degrees

    rng = random.getstate()
    assert all(probability >= 0 for probability in robust_distribution(50, .01, .9))
    get_degrees(50, 100, 11)
    assert random.getstate() == rng
    with pytest.raises(ValueError, match="failure probability"):
        robust_distribution(50, .01, 0)

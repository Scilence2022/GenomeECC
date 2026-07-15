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

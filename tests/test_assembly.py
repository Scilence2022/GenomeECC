"""Forward/reverse assembly with an exact counter and the optional khmer backend."""

from collections import Counter
from random import Random

import pytest

from deBruijnGraph import DeBruijnGraph, assemble_candidates
from ecc_archive import encode_archive, parse_archive, recover_blocks
from fasta import FastaRecord
from utils import rev_seq


def fixture_data():
    sequence = ''.join(Random(908).choices('ACGT', k=480))
    archive = parse_archive(encode_archive(
        [FastaRecord("test", "test topology=linear", sequence, "linear")],
        block_len=80, chunk_size=3, group_blocks=3, rate=.5,
    ))
    return sequence, archive


def exact_counter(segments, k):
    kmers = Counter()
    for segment in segments:
        for start in range(len(segment) - k + 1):
            kmer = segment[start:start + k]
            kmers[min(kmer, rev_seq(kmer))] += 12
    return kmers


def test_crc8_reverse_assembly_crosses_group_boundary():
    sequence, archive = fixture_data()
    graph = DeBruijnGraph(21, counter=exact_counter([sequence[:80], sequence[160:]], 21))
    assert graph.get_kmer_cov(sequence[:21]) == 12
    candidates = assemble_candidates(archive, graph)
    assert set(candidates["test"]) == {0, 2, 3, 4, 5}
    result, = recover_blocks(archive, candidates)
    assert result.sequence == sequence
    assert result.missing_blocks == 1


def test_reverse_assembly_recovers_all_remaining_blocks():
    sequence, archive = fixture_data()
    # No forward anchor coverage in the first chunk. The next chunk's anchor
    # can still reach two blocks backwards before it encounters the first gap.
    graph = DeBruijnGraph(21, counter=exact_counter([sequence[80:]], 21))
    candidates = assemble_candidates(archive, graph)
    assert set(candidates["test"]) == {1, 2, 3, 4, 5}
    assert recover_blocks(archive, candidates)[0].sequence == sequence


def test_partial_terminal_block_uses_crc_of_padded_bases():
    sequence = ''.join(Random(12).choices('ACGT', k=203))
    archive = parse_archive(encode_archive(
        [FastaRecord("tail", "tail", sequence)], block_len=80, chunk_size=3, mode="crc16",
    ))
    graph = DeBruijnGraph(21, counter=exact_counter([sequence], 21))
    assert assemble_candidates(archive, graph)["tail"][2] == sequence[160:]


def test_khmer_fastq_recovery(tmp_path):
    khmer = pytest.importorskip("khmer")
    sequence, archive = fixture_data()
    reads = tmp_path / "synthetic.fq"
    reads.write_text(''.join(
        f"@read-{j}-{i}\n{fragment}\n+\n{'I' * len(fragment)}\n"
        for j, fragment in enumerate([sequence[:80], sequence[160:]]) for i in range(12)
    ))
    graph = DeBruijnGraph(21, counter=khmer.Counttable(21, 1e5, 4))
    graph.count_file(str(reads))
    assert graph.get_kmer_cov(sequence[:21]) == 12
    restored, = recover_blocks(archive, assemble_candidates(archive, graph))
    assert restored.sequence == sequence
    assert restored.missing_blocks == 1

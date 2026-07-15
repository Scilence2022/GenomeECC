"""Multi-chromosome FASTA parsing and genome-level ECC packing."""
import random

import pytest

from genomeECC import GenomeECC


def random_dna(n, seed=0):
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


@pytest.fixture
def genome_fasta(tmp_path):
    """Two NCBI-style records, wrapped at 70 cols, separated by a blank line."""
    c1 = random_dna(4000, 1)
    c2 = random_dna(3000, 2)
    path = tmp_path / "genome.fna"
    with open(path, 'w') as f:
        f.write(">NC_000913.3 Escherichia coli K-12, complete genome circular\n")
        for i in range(0, len(c1), 70):
            f.write(c1[i:i + 70] + "\n")
        f.write("\n")          # blank line: used to truncate the rest of the file
        f.write(">NC_002695.2 Plasmid pX, complete sequence\n")
        for i in range(0, len(c2), 70):
            f.write(c2[i:i + 70] + "\n")
    return path, c1, c2


def test_reads_every_record(genome_fasta):
    path, c1, c2 = genome_fasta
    g = GenomeECC(block_len=320, chunk_size=10)
    g.read_genome(str(path))
    # Standard NCBI headers have no |...| field; extract_id returned None for all of
    # them, collapsing every chromosome onto a single key.
    assert g.chr_num == 2
    assert set(g.genome_seqs) == {'NC_000913.3', 'NC_002695.2'}


def test_sequences_survive_wrapping_and_blank_lines(genome_fasta):
    path, c1, c2 = genome_fasta
    g = GenomeECC(block_len=320, chunk_size=10)
    g.read_genome(str(path))
    assert g.genome_seqs['NC_000913.3'] == c1
    assert g.genome_seqs['NC_002695.2'] == c2


def test_circular_detected_from_description(genome_fasta):
    path, _, _ = genome_fasta
    g = GenomeECC(block_len=320, chunk_size=10)
    g.read_genome(str(path))
    g.genChromesomeECCs()
    assert g.chr_ECC_objs['NC_000913.3'].circular_DNA is True
    assert g.chr_ECC_objs['NC_002695.2'].circular_DNA is False


def test_pipe_delimited_accession_still_wins(tmp_path):
    path = tmp_path / "old.fna"
    path.write_text(">gi|12345|ref|NC_XXXX.1| some description\n" + random_dna(400, 3) + "\n")
    g = GenomeECC(block_len=320, chunk_size=10)
    g.read_genome(str(path))
    assert '12345' in g.genome_seqs


def test_chr_num_available_after_add_chromosome():
    # chr_num was only assigned by read_genome(), so this path raised AttributeError.
    g = GenomeECC(block_len=320, chunk_size=10)
    g.add_chromosome("chr1", random_dna(2000, 4))
    assert g.chr_num == 1
    g.genChromesomeECCs()
    assert len(g.get_all_genome_ecc()) > 0


def test_genome_ecc_is_parseable_back_into_per_chromosome_blobs(genome_fasta):
    path, _, _ = genome_fasta
    g = GenomeECC(block_len=320, chunk_size=10)
    g.read_genome(str(path))
    g.genChromesomeECCs()
    blob = g.get_all_genome_ecc()

    g2 = GenomeECC(block_len=320, chunk_size=10)
    g2.genomeEcc_to_chromosomeECC_arr(blob)
    assert len(g2.chr_ecc_arr) == 2
    for i, head in enumerate(g.chr_ECC_objs):
        assert g2.chr_ecc_arr[i] == g.chr_ECC_objs[head].get_all_ecc_bytes()


def test_ecc_file_round_trip(tmp_path, genome_fasta):
    path, _, _ = genome_fasta
    g = GenomeECC(block_len=320, chunk_size=10)
    g.read_genome(str(path))
    g.genChromesomeECCs()
    out = tmp_path / "genome.ecc"
    g.write_ecc_to_file(str(out))
    assert g.read_ecc_file(str(out)) == g.get_all_genome_ecc()

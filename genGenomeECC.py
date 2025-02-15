from utils import *
from test_utils import *
import zlib
import array
from chromosomeECC import ChromosomeECC
from chromosomeChunk import ChromosomeChunk
from deBruijnGraph import DeBruijnGraph


# seqs = read_fasta('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCF_025643435.1_ASM2564343v1_genomic.fna')


kmer_len = 24
# genome_file = r'genomes/cgl-ssbl-from_nxm.fna'
crc_mode = 'crc16'
# genome_file = r'/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCF_025643435.1_ASM2564343v1_genomic.fna'
#
# genome_file = r'/home/songlf/Git_Repositories/GenomePreservation/genomes/GCA_002847405.1_ASM284740v1_genomic.fna'
# genome_file = r'/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCA_020881995_leaf_2.fa'
genome_file = r'/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/s288c/chr01.fsa'
# genome_file = r'/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/cg13032/cgl-ssbl-from_nxm.fna'

# gc.write_ecc_to_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCF_025643435.1_ASM2564343v1_genomic.fna.ecc')
# gc.write_ecc_to_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCF_025643435.1_ASM2564343v1_genomic.fna.ecc')
# gc.write_ecc_to_file(genome_file + '.ecc')

dbg = DeBruijnGraph()
dbg.set_kmer_len(kmer_len)
print("Counting k-mers file1.....")
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/E0.01-S1-N2E6.r0.read1.fq')
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/E0.01-S1-N2E6.r0.read2.fq')
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/13032_L6_I379.R1.clean.fastq.gz')
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCA_020881995_leaf_2.E0.005-S1-N4E7.r0.read1.fq')
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCA_020881995_leaf_2.E0.005-S1-N4E7.r0.read2.fq')

dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/s288c/chr01.fsa.wgsimc.s.1.fq')
dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/s288c/chr01.fsa.wgsimc.s.2.fq')

# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/s288c/chr01-16-rmt.fsa.wgsimc.1E8.1.fq')
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/s288c/chr01-16-rmt.fsa.wgsimc.1E8.2.fq')

# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/cg13032/cgl-ssbl-from_nxm.fna.wgsimc.s.1.fq')
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/cg13032/cgl-ssbl-from_nxm.fna.wgsimc.s.2.fq')
#
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/s288c/chr01-16-rmt.fsa.fix.reads1.fq')
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/s288c/chr01-16-rmt.fsa.fix.reads2.fq')
# # print("Counting k-mers file2.....")
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/s288c/chr01-16-rmt.fsa.fix.reads2.fq')






dbg.min_kmer_cov = 5
dbg.max_path_num = 1000
dbg.ratio_tolerance = 10
# dbg.block_len = block_len






# block_len = 120
# chunk_size = 5
# redundant_block_rate = 0.04
# kmer_len = 27
# genome_file = r'genomes/cgl-ssbl-from_nxm.fna'
# genome_file = r'/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCF_025643435.1_ASM2564343v1_genomic.fna'
# E0.01-S1-N2E6.r0.read1.fq
#
# genome_file = r'/home/songlf/Git_Repositories/GenomePreservation/genomes/GCA_002847405.1_ASM284740v1_genomic.fna'
# genome_file = r'/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCA_020881995_leaf_2.fa'
# genome_file = r'/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/s288c/chr01-16-rmt.fsa'

# seqs = read_fasta(genome_file)


# dbg.ratio_tolerance = 20
# print("Generating genome ECC codes.....")
# block_len = 200
# chunk_size = 10
# gc = ChromosomeECC(block_len, chunk_size, redundant_block_rate)
# gc.chunk_anchor_byte_num = 8
# # gc.chunk_size = chunk_size
# # gc.set_block_len(block_len)
# gc.set_block_ecc_mode('crc8')
# gc.set_chr_seq(seqs[0])
# gc.get_chunk_crc(0)
#
#
# ecc_len = len(b''.join(gc.chr_ecc()))
# ecc_rate = ecc_len/len(gc.DNA_seq_bytes)

block_len = 320
chunk_size = 10
redundant_block_rate = 0.04

seqs = read_fasta(genome_file)
print("Generating genome ECC codes.....")
gc = ChromosomeECC(block_len, chunk_size, redundant_block_rate)
gc.chunk_anchor_byte_num = 8
# gc.chunk_size = chunk_size
# gc.set_block_len(block_len)
gc.set_block_ecc_mode(crc_mode)
gc.set_chr_seq(seqs[0])
gc.get_chunk_crc(0)
ecc_len = len(b''.join(gc.chr_ecc()))
ecc_rate = ecc_len/len(gc.DNA_seq_bytes)

dbg.ratio_tolerance = 3
dbg.min_kmer_cov = 15
print("Assembling genome blocks .....")
assembled_block_num = 0
failed_blocks = 0


# block_len = 120
# chunk_size = 5
# redundant_block_rate = 0.04
# kmer_len = 24

restored_block_seqs = []

i = 0
failed_block_id = []
too_many_path_num = 0
for i in range(0, gc.chunk_num):
# for i in [100]:
    print('Chunk ID:', end='\t')
    print(i)
    ck = ChromosomeChunk(i)
    ck.set_block_len(gc.block_len)
    ck.chunk_size = gc.chunk_size
    # ck.block_len = block_len
    ck.block_crc_codes = gc.get_chunk_block_crc_codes(i)
    ck.chunk_anchor = gc.get_chunk_anchor(i)
    print(ck.anchor_dna())
    if i + 1 < gc.chunk_num:
        ck.next_chunk_anchor = gc.chunk(i+1)[:6]
    else: # for circular Chromosome
        ck.next_chunk_anchor = gc.chunk(0)[:6]

    block_seqs = dbg.assemble_chunk(ck)
    for j in range(0, len(block_seqs)):
        if block_seqs[j] == gc.chunk_block_seq(i, j):
            assembled_block_num = assembled_block_num + 1
        else:
            failed_block_id.append(i)
            print("Detected issue")
            print(j)
            print(">assembled seq")
            print(block_seqs[j])
            print(">Correct_Block_seq")
            print(gc.chunk_block_seq(i, j))
            failed_blocks = failed_blocks + 1

    print("Block length: " + str(block_len))
    print("Assembled blocks: " + str(assembled_block_num) + " Failed Block num: " + str(failed_blocks))
    print("Block dropout rate: " + str(failed_blocks/(assembled_block_num + failed_blocks)))


for i in range(68, 70):
# for i in [100]:
    print('Chunk ID:', end='\t')
    print(i)
    ck = ChromosomeChunk(i)
    ck.set_block_len(gc.block_len)
    ck.chunk_size = gc.chunk_size
    # ck.block_len = block_len
    ck.block_crc_codes = gc.get_chunk_block_crc_codes(i)
    ck.chunk_anchor = gc.get_chunk_anchor(i)
    print(ck.anchor_dna())
    if i + 1 < gc.chunk_num:
        ck.next_chunk_anchor = gc.chunk(i+1)[:6]
    else: # for circular Chromosome
        ck.next_chunk_anchor = gc.chunk(0)[:6]

    block_seqs = dbg.assemble_chunk(ck)
    for j in range(0, len(block_seqs)):
        if block_seqs[j] == gc.chunk_block_seq(i, j):
            assembled_block_num = assembled_block_num + 1
        else:
            failed_block_id.append(i)
            print("Detected issue")
            print(j)
            print(">assembled seq")
            print(block_seqs[j])
            print(">Correct_Block_seq")
            print(gc.chunk_block_seq(i, j))
            failed_blocks = failed_blocks + 1

    print("Block length: " + str(block_len))
    print("Assembled blocks: " + str(assembled_block_num) + " Failed Block num: " + str(failed_blocks))
    print("Block dropout rate: " + str(failed_blocks/(assembled_block_num + failed_blocks)))

# gc = GenomeECC()
# gc.set_block_len(block_len)
# gc.set_genome_seq(seqs[0])
# gc.get_chunk_crc(0)


dbg1 = DeBruijnGraph()
dbg1.set_kmer_len(kmer_len)
dbg1.count_file(genome_file)
dbg1.count_file(genome_file)
dbg1.count_file(genome_file)
dbg1.min_kmer_cov = 0
dbg1.ratio_tolerance = 20
print("Generating genome ECC codes.....")
block_len = 160
chunk_size = 5
gc = ChromosomeECC(block_len, chunk_size, redundant_block_rate)
gc.chunk_anchor_byte_num = 8
# gc.chunk_size = chunk_size
# gc.set_block_len(block_len)
gc.set_block_ecc_mode('crc8')
gc.set_chr_seq(seqs[5])
gc.get_chunk_crc(0)
ecc_len = len(b''.join(gc.chr_ecc()))
ecc_rate = ecc_len/len(gc.DNA_seq_bytes)
print("Assembling genome blocks .....")
assembled_block_num = 0
failed_blocks = 0

restored_block_seqs = []

i = 0
failed_block_id = []
too_many_path_num = 0
for i in range(0, gc.chunk_num):
# for i in [100]:
    print('Chunk ID:', end='\t')
    print(i)
    ck = ChromosomeChunk(i)
    ck.set_block_len(gc.block_len)
    ck.chunk_size = gc.chunk_size
    # ck.block_len = block_len
    ck.block_crc_codes = gc.get_chunk_block_crc_codes(i)
    ck.chunk_anchor = gc.get_chunk_anchor(i)
    print(ck.anchor_dna())
    if i + 1 < gc.chunk_num:
        ck.next_chunk_anchor = gc.chunk(i+1)[:6]
    else: # for circular Chromosome
        ck.next_chunk_anchor = gc.chunk(0)[:6]

    block_seqs = dbg1.assemble_chunk(ck)
    for j in range(0, len(block_seqs)):
        if block_seqs[j] == gc.chunk_block_seq(i, j):
            assembled_block_num = assembled_block_num + 1
        else:
            failed_block_id.append(i)
            print("Detected issue")
            print(j)
            print(">assembled seq")
            print(block_seqs[j])
            print(">Correct_Block_seq")
            print(gc.chunk_block_seq(i, j))
            failed_blocks = failed_blocks + 1

    print("Block length: " + str(block_len))
    print("Assembled blocks: " + str(assembled_block_num) + " Failed Block num: " + str(failed_blocks))
    print("Block dropout rate: " + str(failed_blocks/(assembled_block_num + failed_blocks)))
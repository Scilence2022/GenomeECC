from utils import *
from test_utils import *
import zlib
import array
from chromosomeECC import ChromosomeECC
from chromosomeChunk import ChromosomeChunk
from deBruijnGraph import DeBruijnGraph

# seqs = read_fasta('cgl-ssbl-from_nxm.fna')
# seqs = read_fasta(r'/home/lifu/Git_Repositories/GenomePreservation/Cg-Genomes/GCA_001447865.2_ASM144786v2_genomic.fna')
# seqs = read_fasta(r'/home/lifu/Git_Repositories/GenomePreservation/Cg-Genomes/GCF_011761195.1_ASM1176119v1_genomic.fna')
# seqs = read_fasta(r'/home/lifu/Git_Repositories/GenomePreservation/Cg-Genomes/' + r'GCF_000404185.1_ASM40418v1_genomic.fna')
# seqs = read_fasta(r'/home/lifu/Git_Repositories/GenomePreservation/Cg-Genomes/' + r'GCA_003025195.1_ASM302519v1_genomic.fna')
#
# seqs = read_fasta(r'/home/lifu/Git_Repositories/GenomePreservation/Cg-Genomes/' + r'GCA_002163915.1_ASM216391v1_genomic.fna')
# seqs = read_fasta(r'/home/lifu/Git_Repositories/GenomePreservation/Cg-Genomes/' + r'GCA_015159635.1_ASM1515963v1_genomic.fna')

# seqs = read_fasta(r'/home/lifu/Git_Repositories/GenomePreservation/Cg-Genomes/' + r'GCF_000196335.1_ASM19633v1_genomic.fna')  # ATCC13032 Wild-type

block_len = 400
chunk_size = 20
redundant_block_rate = 0.04
kmer_cov_cutoff = 4
# genome_file =  r'/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCF_025643435.1_ASM2564343v1_genomic.fna' # r'genomes/cgl-ssbl-from_nxm.fna'

# genome_file = r'/home/songlf/Git_Repositories/GenomePreservation/genomes/GCA_002847405.1_ASM284740v1_genomic.fna' # r'genomes/cgl-ssbl-from_nxm.fna'
genome_file = r'/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/' # r'genomes/cgl-ssbl-from_nxm.fna'


seqs = read_fasta(genome_file)
#__init__(self, block_len=500, chunk_size=6, redundant_block_rate=0.05):

print("Generating genome ECC codes.....")
redundant_block_rate = 0.04
gc = ChromosomeECC(block_len, chunk_size, redundant_block_rate)

# gc.chunk_size = chunk_size
# gc.set_block_len(block_len)
gc.set_chr_seq(seqs[0])
gc.get_chunk_crc(0)
ecc_len = len(b''.join(gc.chr_ecc()))
ecc_rate = ecc_len/len(gc.DNA_seq_bytes)



gc_dec = ChromosomeECC()

unpack = gc_dec.unpack_ecc_bytes(b''.join(gc.chr_ecc()))
gc_dec.reset_with_unpacked_ecc(unpack)
genome_cks = gc_dec.rebuild_chunk_objs_with_ecc_codes(unpack['all_chunk_ecc_bytes'])





dbg = DeBruijnGraph()
print("Counting k-mers.....")
# dbg.count_file(r'/home/lifu/Git_Repositories/GenomePreservation/Cg-Genomes/all/all_cg_genomes.fna')
#dbg.count_file('13032_L6_I379.R1.clean.fastq.cov20')

dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/13032_L6_I379.R1.clean.fastq.cov100')
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/E0.01-S1-N2E6.r0.read1.fq')
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/E0.01-S1-N2E6.r0.read2.fq')

dbg.min_kmer_cov = kmer_cov_cutoff
dbg.max_path_num = 1000
dbg.ratio_tolerance = 20
# dbg.block_len = block_len

print("Assembling genome blocks .....")
assembled_block_num = 0

failed_blocks = 0
restored_block_seqs = []

i = 0
for ck in genome_cks:

    block_seqs = dbg.assemble_chunk(ck)
    for j in range(0, len(block_seqs)):
        if block_seqs[j] == gc.chunk_block_seq(i, j):
            assembled_block_num = assembled_block_num + 1
        else:
            failed_blocks = failed_blocks + 1

    print("Block length: " + str(block_len))
    print("Assembled blocks: " + str(assembled_block_num) + " Failed Block num: " + str(failed_blocks))
    print("Block dropout rate: " + str(failed_blocks/(assembled_block_num + failed_blocks)))
    restored_block_seqs.extend(block_seqs)
    i = i + 1


restored_block_seqs = restored_block_seqs[:gc_dec.block_num]
#Checking dropout blocks
erasure_block_positions = []
for i in range(0, len(restored_block_seqs)):
    if restored_block_seqs[i] == "":
        restored_block_seqs[i] = "A" * block_len
        erasure_block_positions.append(i)


gc_dec.block_seqs = restored_block_seqs
print("Restoring block bytes ..... ")
gc_dec.block_seqs_to_bytes()
gc_dec.block_bytes_to_colum_bytes()
print("Restoring the missing blocks by the RS codes ..... ")
gc_dec.restore_colum_bytes_rs(unpack['all_rs_code_bytes'], erasure_block_positions)
gc_dec.colum_bytes_to_block_bytes()
succ_num = 0
for i in range(0, gc_dec.block_num):
    if gc.block_bytes[i] == gc_dec.block_bytes[i]:
        succ_num = succ_num + 1
    else:
        print(i, end='\t')

assembled_genome = gc_dec.assembled_chr_seq()

if assembled_genome == seqs[0]:
    print("Genome decoded correctly without any mistakes!")
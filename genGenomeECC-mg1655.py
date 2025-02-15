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
genome_file = r'/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCF_025643435.1_ASM2564343v1_genomic.fna'
genome_file = r'/data/songlf/2.GenomeCali/05.correct_song_result_20250210/03.WRX_MG1655/R.MG1655.fasta'
genome_file = r'/data/songlf/2.GenomeCali/05.correct_song_result_20250210/03.WRX_MG1655/01.result.WRX_MG1655.fasta'
genome_file = r'/data/songlf/2.GenomeCali/大肠谷棒酵母基因组数据/R.MG1655-Fixed.fa'


#


dbg = DeBruijnGraph()
dbg.set_kmer_len(kmer_len)
print("Counting k-mers file1.....")
# dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCF_025643435.1_ASM2564343v1_genomic.fna.wgsimc.1E6.1.fq')
# dbg.count_file('/data/songlf/2.GenomeCali/大肠谷棒酵母基因组数据/WRX_MG1655_S4_R1_001.fastq')
# print("Counting k-mers file2.....")
# # dbg.count_file('/data/songlf/0.DNA_Storage/0.GenomePreservation/genomes/GCF_025643435.1_ASM2564343v1_genomic.fna.wgsimc.1E6.2.fq')
# dbg.count_file('/data/songlf/2.GenomeCali/大肠谷棒酵母基因组数据/WRX_MG1655_S4_R2_001.fastq')





# dbg.block_len = block_len


block_len = 320
chunk_size = 10
redundant_block_rate = 0.06

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





dbg.min_kmer_cov = 2
dbg.ratio_tolerance = 20
dbg.max_path_num = 1000

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
# for i in [1234]:
# for i in range(0, 100):
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
            print("Checking blocks: " + str(j))
            print(test_strand_in_dbg(gc.chunk_block_seq(i, j), dbg))

            # print("Detected issue")
            # print(j)
            # print(">assembled seq")
            # print(block_seqs[j])
            # print(">Correct_Block_seq")
            # print(gc.chunk_block_seq(i, j))
            failed_blocks = failed_blocks + 1

    print("Block length: " + str(block_len))
    print("Assembled blocks: " + str(assembled_block_num) + " Failed Block num: " + str(failed_blocks))
    print("Block dropout rate: " + str(failed_blocks/(assembled_block_num + failed_blocks)))


gc_dec = ChromosomeECC()

unpack = gc_dec.unpack_ecc_bytes(b''.join(gc.chr_ecc()))
gc_dec.reset_with_unpacked_ecc(unpack)
genome_cks = gc_dec.rebuild_chunk_objs_with_ecc_codes(unpack['all_chunk_ecc_bytes'])



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




#
# for i in range(68, 70):
# # for i in [100]:
#     print('Chunk ID:', end='\t')
#     print(i)
#     ck = ChromosomeChunk(i)
#     ck.set_block_len(gc.block_len)
#     ck.chunk_size = gc.chunk_size
#     # ck.block_len = block_len
#     ck.block_crc_codes = gc.get_chunk_block_crc_codes(i)
#     ck.chunk_anchor = gc.get_chunk_anchor(i)
#     print(ck.anchor_dna())
#     if i + 1 < gc.chunk_num:
#         ck.next_chunk_anchor = gc.chunk(i+1)[:6]
#     else: # for circular Chromosome
#         ck.next_chunk_anchor = gc.chunk(0)[:6]
#
#     block_seqs = dbg.assemble_chunk(ck)
#     for j in range(0, len(block_seqs)):
#         if block_seqs[j] == gc.chunk_block_seq(i, j):
#             assembled_block_num = assembled_block_num + 1
#         else:
#             failed_block_id.append(i)
#             print("Detected issue")
#             print(j)
#             print(">assembled seq")
#             print(block_seqs[j])
#             print(">Correct_Block_seq")
#             print(gc.chunk_block_seq(i, j))
#             failed_blocks = failed_blocks + 1
#
#     print("Block length: " + str(block_len))
#     print("Assembled blocks: " + str(assembled_block_num) + " Failed Block num: " + str(failed_blocks))
#     print("Block dropout rate: " + str(failed_blocks/(assembled_block_num + failed_blocks)))
#
# # gc = GenomeECC()
# # gc.set_block_len(block_len)
# # gc.set_genome_seq(seqs[0])
# # gc.get_chunk_crc(0)
#
#
# dbg1 = DeBruijnGraph()
# dbg1.set_kmer_len(kmer_len)
# dbg1.count_file(genome_file)
# dbg1.count_file(genome_file)
# dbg1.count_file(genome_file)
# dbg1.min_kmer_cov = 0
# dbg1.ratio_tolerance = 20
# print("Generating genome ECC codes.....")
# block_len = 160
# chunk_size = 5
# gc = ChromosomeECC(block_len, chunk_size, redundant_block_rate)
# gc.chunk_anchor_byte_num = 8
# # gc.chunk_size = chunk_size
# # gc.set_block_len(block_len)
# gc.set_block_ecc_mode('crc8')
# gc.set_chr_seq(seqs[5])
# gc.get_chunk_crc(0)
# ecc_len = len(b''.join(gc.chr_ecc()))
# ecc_rate = ecc_len/len(gc.DNA_seq_bytes)
# print("Assembling genome blocks .....")
# assembled_block_num = 0
# failed_blocks = 0
#
# restored_block_seqs = []
#
# i = 0
# failed_block_id = []
# too_many_path_num = 0
# for i in range(0, gc.chunk_num):
# # for i in [100]:
#     print('Chunk ID:', end='\t')
#     print(i)
#     ck = ChromosomeChunk(i)
#     ck.set_block_len(gc.block_len)
#     ck.chunk_size = gc.chunk_size
#     # ck.block_len = block_len
#     ck.block_crc_codes = gc.get_chunk_block_crc_codes(i)
#     ck.chunk_anchor = gc.get_chunk_anchor(i)
#     print(ck.anchor_dna())
#     if i + 1 < gc.chunk_num:
#         ck.next_chunk_anchor = gc.chunk(i+1)[:6]
#     else: # for circular Chromosome
#         ck.next_chunk_anchor = gc.chunk(0)[:6]
#
#     block_seqs = dbg1.assemble_chunk(ck)
#     for j in range(0, len(block_seqs)):
#         if block_seqs[j] == gc.chunk_block_seq(i, j):
#             assembled_block_num = assembled_block_num + 1
#         else:
#             failed_block_id.append(i)
#             print("Detected issue")
#             print(j)
#             print(">assembled seq")
#             print(block_seqs[j])
#             print(">Correct_Block_seq")
#             print(gc.chunk_block_seq(i, j))
#             failed_blocks = failed_blocks + 1
#
#     print("Block length: " + str(block_len))
#     print("Assembled blocks: " + str(assembled_block_num) + " Failed Block num: " + str(failed_blocks))
#     print("Block dropout rate: " + str(failed_blocks/(assembled_block_num + failed_blocks)))
import random
from utils import *
from DNAdroplet import DNADroplet
#from fountain import Fountain
from DNAfountain import DNAFountain
from glass import Glass
from crc16pure import *
from test_utils import *
import copy

from reedsolo import RSCodec, ReedSolomonError
# rsc = RSCodec(10)  # 10 ecc symbols


delta = 0.01
c_value = 0.01

data_block_length = 30
test_num = 100
fountain_seed = 2
p1 = 'CCTGCAGAGTAGCATGTC'  # 5'-->3'
p2 = 'CTGACACTGATGCATCCG'  # complement seq of P2


work_dir = r'.'  + '/'

input_file = work_dir + r'cgl-ssbl-from_nxm.fna'
afile = open(input_file, 'rb')

out_handle = "future-10-jpg_seed2_t100"

genome_seqs = read_fasta(input_file)
genome_seqs_bytes = DNA_seqs_to_bytes(genome_seqs)


filebytes = b''


passwd = 0b10101010101111011001010100110110

fountain_init_index = 211011011
ft = DNAFountain(filebytes, data_block_length, fountain_init_index, fountain_seed, 16, 16)
ft.fix_bytes()
ft.des = True
ft.sec_key = passwd.to_bytes(8, byteorder ='big')

total_size = 12000 #int(fdna1.num_of_chunks * 1.15)
core_size = int(total_size*0.98)

# droplet_four_pics = []

droplet_all = get_droplets(total_size, ft)

chunk_size = ft.chunk_size
chunk_num = ft.num_of_chunks
degree_table = get_degrees(chunk_num, chunk_num * 3, fountain_seed, delta, c_value)
# cup = Glass(chunk_num)

# j = 1
# adp = DNADroplet(bytes(chunk_size))
# adp.des = True
# adp.sec_key = ft.sec_key
# adp.num_of_chunks = chunk_num
# adp.set_index(droplet_all[j].index)
# print(adp.set_droplet_from_DNA(droplet_all[j].to_DNA()))
# adp1 = droplet_all[j]



file2 = open(input_file + r'.sim.'  + out_handle, 'tw')
#file2.write('Head Index\tData\tDNA\tDNA-Primers\tDegree\tChunk Nums\tTail Index\n')
for dps in droplet_all:
    file2.write(str(dps.index))
    file2.write('\t')
    file2.write(p1)
    file2.write(dps.to_DNA())
    file2.write(p2)
    file2.write('\n')
file2.close()



file2 = open(input_file + r'.tab.rich.' + out_handle, 'tw')
file2.write('Head Index\tData\tDNA\tDNA-Primers\tDegree\tChunk Nums\tTail Index\n')
for dps in droplet_all:
    file2.write(str(dps.index))
    file2.write('\t')
    file2.write(str(dps.data))
    file2.write('\t')
    file2.write(dps.to_DNA())
    file2.write('\t')
    file2.write(p1)
    file2.write(dps.to_DNA())
    file2.write(p2)
    file2.write('\t')
    file2.write(str(dps.degree))
    file2.write('\t')
    file2.write(str(dps.get_chunk_nums()))
    file2.write('\n')
file2.close()


# droplet_four_pics.append(droplet_all)
suc_num = 0
for i in range(0, test_num):
    droplet_sample = []
    # droplet_sample[i] = random.sample(range(0, total_size), core_size)
    random.seed(i + 1)

    for j in random.sample(range(0, total_size), core_size):
        adp = DNADroplet(bytes(chunk_size))
        adp.des = True
        adp.sec_key = ft.sec_key
        adp.num_of_chunks = chunk_num
        adp.set_index(droplet_all[j].index)
        adp.set_droplet_from_DNA(droplet_all[j].to_DNA())
        adp.degree = degree_table[droplet_all[j].index % len(degree_table)]
        adp.update()
        droplet_sample.append(adp)
        #droplet_sample.append(copy.deepcopy(droplet_all[j]))


    if test_droplets(droplet_sample, ft):
        suc_num = suc_num + 1

    print(suc_num, end="/")
    print(i + 1)
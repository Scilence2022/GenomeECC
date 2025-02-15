import sys; print('Python %s on %s' % (sys.version, sys.platform))

from utils import *
from test_utils import *
from deBruijnGraph import DeBruijnGraph
from glass import Glass
from DNAdroplet import DNADroplet
import time
import os.path
import getopt

delta = 0.01
c_value = 0.01

input_file = r'input_files/future-10.jpg.sim.future-10-jpg_seed2_t100.fasta'
# file_type = 'dump_kmers'
output_file = 'output-10.jpg'
kmer_size = 31
kmer_cut_off = 0
chunk_size = 30
chunk_num = 10517
fountain_seed = 2
index_bytes = 4
crc_bytes = 4
# double_index = False
# both_way = True
index_l = 211011011
index_u = (chunk_num * 2) + index_l
sec_key = 0b10101010101111011001010100110110
# sec_key = 0b10101010101111011001010100110111 # wrong key
sec_key = sec_key.to_bytes(8, byteorder ='big')

def droplet_tp():
    adp = DNADroplet(bytes(chunk_size))
    adp.des = True
    adp.sec_key = sec_key
    adp.num_of_chunks = chunk_num
    return adp




opts,args = getopt.getopt(sys.argv[1:],'-h-i:-o:t:-k:-c:-n:-s:-a-b',
                          ['help','input=','output=', 'file_type=', 'kmer_size=', 'chunk_size=', 'chunk_num=', 'seed=', 'anchor', 'both_way', 'min_index=',  'max_index=', 'index_bytes=', 'ec_bytes='])

usage = 'Usage:\n' + r'      python decode.py -i input_file -t type_of_seqs -o outfile [Options]'
options = 'Options:\n'
options = options + r'      -h, --help                              Show help information' + '\n'
options = options + r'      -i, --input   <input file>              Input file' + '\n'
options = options + r'      -t, --file_type   <file type>           Input file type: FastQ, Fasta or Jellyfish dumped k-mers (default)' + '\n'
options = options + r'      -o, --output  <output file>             Output file' + '\n'
options = options + r'      -k, --kmer_size  <number>               k-mer size, default = 21 ' + '\n'
options = options + r'      -c, --chunk_size  <size>                Chunk size, default = 30 (bytes)' + '\n'
options = options + r'      -n, --chunk_num  <number>               Chunk number, default = 10,741 (for testing only)' + '\n'
options = options + r'      --cut  <number>                         Cut_off for elimination of low coverage k-mers, default = 0 ' + '\n'
options = options + r'      -s, --seed    <seed>                    Fountain random seed, default 1' + '\n'
# options = options + r'      -a, --anchor                            Anchor codes, default On' + '\n'
# options = options + r'      -b, --both_way                          Both-way search mode, default On' + '\n'
options = options + r'      --min_index  <initial index>            Initial index, default = 1' + '\n'
options = options + r'      --max_index  <max index>                Max index, default = 20000' + '\n'
options = options + r'      --index_bytes  <number>                 Length of index and anchor codes, default = 4 (bytes)' + '\n'
options = options + r'      --ec_bytes  <number>                    Length of ec codes, default = 4 (bytes)' + '\n'


for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print(usage)
        print(options)
        sys.exit()
    if opt_name in ('-i','--input'):
        input_file = opt_value
    if opt_name in ('-o','--output'):
        output_file = opt_value
    if opt_name in ('-t','--file_type'):
        file_type = opt_value
    if opt_name in ('-k','--kmer_size'):
        kmer_size = int(opt_value)
    if opt_name in ('-c','--chunk_size'):
        chunk_size = int(opt_value)
    if opt_name in ('-n', '--chunk_num'):
        chunk_num = int(opt_value)
    if opt_name in ('-s', '--seed'):
        f_seed = int(opt_value)
    # if opt_name in ('-a', '--anchor'):
    #     double_index = True
    if opt_name in ('--min_index', '--notmatch'):
        index_l = int(opt_value)

    if opt_name in ('--max_index', '--notmatch'):
        index_u = int(opt_value)

    if opt_name in ('--index_bytes', '--notmatch'):
        index_bytes = int(opt_value)

    if opt_name in ('--crc_bytes', '--notmatch'):
        crc_bytes = int(opt_value)


start = time.perf_counter()

deG = DeBruijnGraph()
deG.set_kmer_len(kmer_size)

print('\nCounting k-mers ......')
deG.count_file(input_file)

# deG.kmer_len = kmer_size
# deG.veri_kmers = False
# deG.max_path_num = 10000

# if file_type == 'FastQ' or file_type == 'fastq':
#     deG.open_fastq(input_file)
# else:
#     if file_type == 'fasta' or file_type == 'Fasta':
#         deG.open_fasta(input_file)
#     else:
#         deG.open_dump(input_file)

print('\nReconstructing DNA droplets ......')

degree_table = get_degrees(chunk_num, chunk_num * 3, fountain_seed, delta, c_value)

a = time.perf_counter()
cup = Glass(chunk_num)

droplet_template = droplet_tp()
for index in range(index_l, index_u+1):
    fd_droplets = deG.find_droplets(index, droplet_template)
    if len(fd_droplets) == 1:
        fd_droplets[0].degree = degree_table[index % len(degree_table)]
        fd_droplets[0].update()
        cup.addDroplet(fd_droplets[0])
    else:
        if len(fd_droplets) > 1:
            print("Multiple droplets found with index " + str(index))



print('Decoding by fountain codes .........')

cup.decode()
if cup.isDone():
    cup.writeToFile(output_file)
else:
    print("Failed to decode the original information.")

end = time.perf_counter() - start
print('Decoding time: ', end ='')
print(end)




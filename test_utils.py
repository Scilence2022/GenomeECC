import crc16pure
# import random
# import fountain
# import droplet
import copy
from DNAdroplet import DNADroplet
from utils import *
from glass import Glass
from DNAfountain import DNAFountain
import re
import operator
import math
import glass

#20190701
def get_droplets(num, fdna1):
    i = 0
    droplets = []
    gc_drop_num = 0
    adrop = None
    print(i)
    print(num)
    while i < num:
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA()):
            droplets.append(adrop)
            # print(str(i)+"     Good Droplet detected!")
            i = i + 1
        else:
            gc_drop_num = gc_drop_num + 1
        #     print("Bad Droplet detected!")
        #     print(i)
    print('GC drop num:',end='\t')
    print(gc_drop_num)
    return droplets

def get_droplets_check_repeat_kmer(num, fdna1, kmer_len=21):
    i = 0
    droplets = []
    #kmer_lenth = kmer_len - 1
    kms = {}
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):
            # if deGraph.highst_km_freq(adrop.to_DNA_CRC()) < 5:
            dps_kms = kmers_of_str(adrop.to_DNA_CRC(), kmer_len - 1)
            # kmers_in_dict()
            if not any_kmers_in_dict(dps_kms, kms):
                droplets.append(adrop)
                i = i + 1
                # Adding droplet kmers into kms
                for km in dps_kms:
                    kms[km] = 1
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

    #     if j > 10000:
    #         j = 0
    #         print(gc_drop_num, end='\t')
    #         print(km_drop_num, end='\t')
    #         print(num)
    # print('GC drop num:', end='\t')
    # print(gc_drop_num, end='\t')
    # print('Km drop num:', end='\t')
    # print(km_drop_num, end='\t')
    # print(num)
    return droplets



def get_droplets_check_repeat_kmer_multi_ft(num, fdna1, deGraph):

    kmer_len = deGraph.kmer_len
    i = 0
    droplets = []
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):

            dps_kms = kmers_of_str(adrop.to_DNA_CRC(), kmer_len - 1)

            if not any_kmers_in_dict(dps_kms, deGraph.kmers):
                droplets.append(adrop)
                i = i + 1
                for km in dps_kms:
                    deGraph.kmers[km] = 1
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

        if j > 10000:
            j = 0
            print(len(droplets), end='\t')
            print(gc_drop_num, end='\t')
            print(km_drop_num, end='\t')
            print(num)
    print(len(droplets), end='\t')
    print(gc_drop_num, end='\t')
    print(km_drop_num, end='\t')
    print(num)
    return droplets


def get_droplets_check_repeat_kmer_4deG(num, fdna1, deGraph, max_kmer_repeat_num=5):

    kmer_len = deGraph.kmer_len
    i = 0
    droplets = []
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):
            if deGraph.highst_km_freq(adrop.to_DNA_CRC()) <= max_kmer_repeat_num:
                droplets.append(adrop)
                i = i + 1
                deGraph.add_seq(adrop.to_DNA_CRC())
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

        if j > 10000:
            j = 0
            print(gc_drop_num, end='\t')
            print(km_drop_num, end='\t')
            print(num)
    print(gc_drop_num, end='\t')
    print(km_drop_num, end='\t')
    print(num)
    return droplets

def add_droplets_to_deG_dist(num, fdna1, deGraph, min_kmer_dist=1):
    # kmer_len = deGraph.kmer_len
    i = 0
    droplets = []
    gc_drop_num = 0
    km_drop_num = 0
    # adrop = None
    j = 0
    while i < num:
        j = j + 1
        adrop = fdna1.DNAdroplet()
        if check_dna(adrop.to_DNA_CRC()):
            if kmer_editing_dist_seq_deG(adrop.to_DNA_CRC(), deGraph) >= min_kmer_dist:
                droplets.append(adrop)
                i = i + 1
                deGraph.add_seq(adrop.to_DNA_CRC())
            else:
                km_drop_num = km_drop_num + 1
        else:
            gc_drop_num = gc_drop_num + 1

        if j > 10000:
            j = 0
            print(gc_drop_num, end='\t')
            print(km_drop_num, end='\t')
            print(num)
    print(gc_drop_num, end='\t')
    print(km_drop_num, end='\t')
    print(num)
    return droplets


def dropout_rate(num, fdna1, min_gc=0.45, max_gc=0.55, max_homo=5, max_repeat_kmer=5):
    i = 0
    # droplets = []
    # deGraph = DeBruijnGraph()
    drop_result = {}
    drop_result['gc'] = 0
    drop_result['homo'] = 0

    while i < num:
        adrop = fdna1.DNAdroplet()
        dnastr = adrop.to_DNA_CRC()
        gc_rate = calc_gc(dnastr)
        if gc_rate > max_gc:
            drop_result['gc'] = drop_result['gc'] + 1
        if gc_rate < min_gc:
            drop_result['gc'] = drop_result['gc'] + 1
        homo_poly_len = max_homo_len(dnastr)
        if homo_poly_len > max_homo:
            drop_result['homo'] = drop_result['homo'] + 1
#         if check_dna(adrop.to_DNA_CRC()):
#             # print(str(i)+"     Good Droplet detected!")
# #           print(deGraph.highstKmFreq(adrop.to_DNA_CRC()))
#             droplets.append(adrop)
#             # deGraph.addSeq(adrop.to_DNA_CRC())
#             #print("Good Droplet detected!")
#
#         else:
#             drop_num = drop_num + 1
        i = i + 1
    return drop_result

def get_degree_droplet(degree, fdna1):
    adrop = fdna1.DNAdroplet()
    while adrop.degree != degree:
        adrop = fdna1.DNAdroplet()

    return adrop

def test_glass(gt, fdna1):
    i = 0
    while i < gt.num_chunks - 1:
        if gt.chunks[i] != fdna1.chunk(i):
            # print(i)
            return False
        i = i + 1
    gap = len(gt.chunks[i]) - len(fdna1.chunk(i))
    if gt.chunks[i] != fdna1.chunk(i) + b'\x00' * gap:
        # print(i)
        return False
    return True


def test_droplets(droplets, fdna1):
    gt = glass.Glass(fdna1.num_of_chunks)
    i = 0
    for droplet in droplets:
        # print(str(i) + "  adding droplets")
        # print(droplet.degree)
        if droplet.des:
            droplet.decry_data()
        gt.addDroplet(droplet)
        i = i + 1
        #print(i)
    gt.decode()
    i = 0
    while i < gt.num_chunks-2:
        if gt.chunks[i] != fdna1.chunk(i):
            # print(i)
            # print(gt.chunks[i])
            # print(fdna1.chunk(i))
            return False
        i = i + 1
    gap = len(gt.chunks[i]) - len(fdna1.chunk(i))
    if gt.chunks[i] != fdna1.chunk(i) + b'\x00' * gap:
        #print(i)
        return False
    return True

def test_DNA_collection(DNAs, fdna1, kmer_len=15 , min_index=1, max_index=100000):

    degree_table = get_degrees(fdna1.num_of_chunks, int(fdna1.num_of_chunks * 10), fdna1.seed)
    deG = DeBruijnGraph()
    deG.kmer_len = kmer_len
    deG.add_seqs(DNAs)
    test_results = {}

    print("Recovering data...")
    i = min_index
    find_path_num = 0
    find_multi_path_num = 0
    full_recover = False
    aGlass = Glass(fdna1.num_of_chunks)
    recov_droplets = []
    while (i <= max_index):
        # print("Finding index " + str(i) + " path")
        deG.find_droplet_DNA(i, 19)
        if (len(deG.crc_paths) > 0):
            find_path_num = find_path_num + 1
            # print("findpathnum  " + str(find_path_num))
            if (len(deG.crc_paths) > 1):
                find_multi_path_num = find_multi_path_num + 1
            else:
                adrop = DNADroplet()
                adrop.num_of_chunks = fdna1.num_of_chunks
                adrop.degree = degree_table[i - 1]
                adrop.set_droplet_from_DNA_CRC(deG.crc_paths[0])
                recov_droplets.append(copy.deepcopy(adrop))
                aGlass.addDroplet(adrop)
        print("Searched index: " + str(i) + " Found Droplets: " + str(find_path_num) + " Found multiple Paths: " + str(find_multi_path_num))
        i = i + 1
    if (test_glass(aGlass, fdna1)):
        full_recover = True

    test_results['Droplet_found'] = find_path_num
    test_results['Multi_path_found'] = find_multi_path_num
    test_results['Full_recover'] = full_recover

    return test_results


def recover_file(droplets, num_of_chunks,output_file,fdna1):
    OUT = open(output_file, 'wb')
    error_chunk_num = 0
    gt = Glass(num_of_chunks)
    i = 0
    for droplet in droplets:
        gt.addDroplet(droplet)
        i = i + 1
        #print(i)
    i = 0
    while i < gt.num_chunks-1:
        OUT.write(gt.chunks[i])
        if gt.chunks[i] != fdna1.chunk(i):
            #print(i)
            error_chunk_num = error_chunk_num + 1
        i = i + 1

    gap = len(gt.chunks[i]) - len(fdna1.chunk(i))
    OUT.write(gt.chunks[i])
    if gt.chunks[i] != fdna1.chunk(i) + b'\x00' * gap:
        error_chunk_num = error_chunk_num + 1

    return error_chunk_num


def get_droplets_noCheck(num, fdna1):
    i = 0
    droplets = []
    adrop = None
    while i < num:
        adrop = fdna1.DNAdroplet()
        #if check_dna(adrop.to_DNA_CRC()):
        droplets.append(adrop)
            #print("Good Droplet detected!")
        i = i + 1
        #else:
            #print("Bad Droplet detected!")
         #print(i)
    return droplets

def count_length_DNAs(dnas, max_len = 162):
    lengths = {}
    for j in range(0, max_len+1):
        lengths[j] = 0
    for dna in dnas:
        lengths[len(dna)] = lengths[len(dna)] + 1
    return lengths

def count_avg_length_DNAs(dnas):
    num_of_DNAs = len(dnas)
    length_all_DNAs = 0
    for dna in dnas:
        length_all_DNAs = length_all_DNAs + len(dna)
    return int(length_all_DNAs/num_of_DNAs)


# 2020/06/30
def hashToFile(hs, file):
    out = open(file, 'tw')
    for a in hs:
        out.write(str(a))
        out.write("\t")
        out.write(str(hs[a]))
        out.write("\n")
    out.close()


# Max value of hash
def max_value_hash(hs):
    max = 0
    for aKey in hs:
        if hs[aKey] > max:
            max = hs[aKey]
    return max

def poisson_seq_num(ranmd, rep_f_exp, rep_times=20, size=1):
     init_seq_nums = np.random.poisson(ranmd)
     rep_folds = math.pow((float(1) + poisson_rep_factor(rep_f_exp)), rep_times)
     return int(init_seq_nums * rep_folds)

def poisson_rep_factor(exp):
    ranmd = int(exp*10)
    rep_f = float(np.random.poisson(ranmd)/10)
    if rep_f > 1:
        rep_f = 1
    return rep_f


def key_kmers_of_strand(dnastr, kmer_len, overlap=20):
    kmers = []
    str_len=len(dnastr)
    kmers.append(kmers_of_position(dnastr, kmer_len, 0))
    kmers.append(kmers_of_position(dnastr,kmer_len,overlap))
    kmers.append(kmers_of_position(dnastr, kmer_len, str_len-overlap-kmer_len))
    kmers.append(kmers_of_position(dnastr, kmer_len, str_len - kmer_len))
    return kmers

    # 20201218 read aln file
def read_fasta(file):
    f = open(file, "r")
    seqs = []
    matchLineA = re.compile('^(>)')

    line = f.readline()

    seq = ''
    while line.strip():
        if not matchLineA.match(line):
            seq = seq + line.strip()
        else:
            if seq != '':
                seqs.append(seq)
                seq = ''
        line = f.readline()
    f.close()

    if seq != '':
        seqs.append(seq)
    return seqs


def hash_keys_values(hs):
    for a in hs:
        print(a, end='\t')
        print(hs[a])


def sub_seqs(seqs, ff=0, ee=-1):
    if ee < 0:
        ee = len(seqs[0]) -1
    arr = []
    for a in seqs:
        arr.append(a[ff:ee])
    return arr


def kmer_editing_dist(k1, k2):
    assert len(k1) == len(k2), "length of kmer 1 and 2 are not same"
    dis = 0
    for i in range(0,len(k1)):
        if k1[i] != k2[i]:
            dis = dis + 1
    return dis


def split_str(str, block_len=1000):
    str_len = len(str)
    block_num = math.ceil(str_len/block_len)
    split_strs = []
    for i in range(0, block_num):
        if i < block_num:
            split_strs.append(str[i*block_len:i*block_len + block_len])
        else:
            split_strs.append(str[i*block_len:])
    return split_strs


def seqs_to_fasta_file(seqs, file):
# 2020/06/30
    out = open(file,'tw')
    for i in range(0, len(seqs)):
        out.write(">" + str(i) + "\n")
        out.write(str(seqs[i]))
        out.write("\n")
    out.close()


def introduce_random_errors_array(arr, error_rate):
    """
    Introduce random errors to an array of integers.

    Parameters:
        arr (list of int): The array to introduce errors to.
        error_rate (float): The rate of errors to introduce, as a decimal value between 0 and 1.

    Returns:
        list of int: The array with errors introduced.
    """
    new_arr = []
    for num in arr:
        if random.random() < error_rate:
            # introduce error
            new_arr.append(random.randint(0, max(arr)))
        else:
            new_arr.append(num)
    return new_arr


def introduce_random_errors_bytes(data: bytes, error_rate: float) -> bytes:
    """
    Introduces random errors to a bytes object.

    Args:
        data (bytes): The bytes object to introduce errors to.
        error_rate (float): The error rate, as a fraction between 0 and 1.

    Returns:
        bytes: A new bytes object with errors introduced.
    """
    num_errors = int(len(data) * error_rate)
    error_indices = random.sample(range(len(data)), num_errors)
    result = bytearray(data)
    for i in error_indices:
        result[i] = random.randint(0, 255)
    return bytes(result)


def int_array_to_bytes(int_array, byte_size = 1) -> bytes:
# def int_array_to_bytes(int_array: list[int], byte_size=1) -> bytes:
    """
    Converts an array of integers to a bytes object.

    Args:
        int_array (List[int]): The array of integers to convert.

    Returns:
        bytes: A bytes object representing the same data as the input array.
    """
    byte_list = [i.to_bytes(byte_size, byteorder='big') for i in int_array]
    return b''.join(byte_list)


def bytes_to_int_array(data, chunk_size):
    """
    Encode a bytes object into an array of integers by splitting the bytes with a specific length.

    Args:
        data (bytes): The bytes object to encode.
        chunk_size (int): The length of each chunk.

    Returns:
        List[int]: An array of integers representing the encoded data.
    """
    encoded_data = []
    for i in range(0, len(data), chunk_size):
        chunk = data[i:i + chunk_size]
        value = int.from_bytes(chunk, byteorder='big')
        encoded_data.append(value)
    return encoded_data


def test_strand_in_dbg(dnastr, dbg):
    km_list = kmer_list_of_str(dnastr,dbg.kmer_len)
    min_cov = dbg.min_kmer_cov
    all_km_num = len(km_list)
    miss_km_num = 0
    miss_km_hs = {}
    for index, item in enumerate(km_list):
        print(item, end='\t')
        print(dbg.get_kmer_cov(item), end ='  ||  ')
        if dbg.get_kmer_cov(item) < min_cov:
            miss_km_num = miss_km_num + 1
            print(item)
    return miss_km_num





# def test_ch_chunk()
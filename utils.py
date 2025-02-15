import random
import math
import numpy as np
from random import choices
from pyDes import des, CBC, PAD_PKCS5, ECB
import sha256

def	charN(str, N):
	if N < len(str):
		return str[N]
	return 'X'


def oxor(str1, str2):
	length = max(len(str1), len(str2))
	return ''.join(chr(ord(charN(str1, i)) ^ ord(charN(str2, i))) for i in range(length))

def	xor(str1, str2):
	length = max(len(str1), len(str2))
	if len(str1) > len(str2):
		str2 = str2 + bytes(len(str1) - len(str2))
	if len(str2) > len(str1):
		str1 = str1 + bytes(len(str2) - len(str1))
	allBytes = b''
	i = 0
	while i < length:
		allBytes = allBytes + bytes([str1[i] ^ str2[i]])
		i =  i + 1
	return allBytes


def randChunkNums(degree, num_chunks):
	size = random.randint(1,max(5, int(num_chunks*3/5)))
	return random.sample(range(num_chunks), size)


def get_degrees(N, k, symbol_index=-1, delta=0.1, c_value=0.2, distribution_name='robust'):
	#print(distribution_name)
	if distribution_name == "ideal":
		probabilities = ideal_distribution(N)
	elif distribution_name == "robust":
		probabilities = robust_distribution(N, c_value, delta)
	else:
		probabilities = None

	population = list(range(0, N))
	if symbol_index > 0:
		random.seed(symbol_index)
	#Not initialized yet
	else:
		return -1
	return choices(population, probabilities, k=k)


def generate_chunk_nums(blocks_quantity, degree, symbol_index):
	random.seed(symbol_index)
	indexes = random.sample(range(blocks_quantity), degree[0])
	return indexes


def ideal_distribution(K):

	probabilities = [0, 10 / K]
	probabilities += [1 / (k * (k - 1)) for k in range(2, K)]
	probabilities_sum = np.sum(probabilities)
	probabilities /= probabilities_sum

	# assert probabilities_sum >= 1 - epsilon and probabilities_sum <= 1 + epsilon, "The ideal distribution should be standardized"
	return probabilities

def robust_distribution(K, c_value = 0.5, robust_failure_probability=0.1):

	# print('c value=', end='\t')
	# print(c_value, end='\t')
	# print('epsilon', end = '\t')
	# print(robust_failure_probability)

	S = c_value * math.log(K / robust_failure_probability) * math.sqrt(K)
	print(S)
	# M = round(K/S)
	M = round(K / S)
	print(M)
	extra_proba = [0] + [1/(i * M) for i in range(1, M-1)]
	extra_proba += [S * math.log(S / robust_failure_probability) / K]  # Spike at M
	extra_proba += [0 for k in range(M, K)]

	probabilities = np.add(extra_proba, ideal_distribution(K))
	print(np.sum(probabilities))
	probabilities /= np.sum(probabilities)
	#probabilities_sum = np.sum(probabilities)
	#assert probabilities_sum >= 1 - epsilonc and probabilities_sum <= 1 + epsilonc, "The robust distribution should be standardized"
	return probabilities

def DNA_seqs_to_bytes(seqs):
	byte_seqs = []
	for sq in seqs:
		byte_seqs.append(DNAToBytes(sq))
	return byte_seqs


def bytesToDNA(manyBytes):
	dnastr = ''
	for aByte in manyBytes:
		dnastr = dnastr + byteToDNA(aByte)
	return dnastr


def DNAToBytes(dnaStr):
	i = 0
	nBytes=b''
	dnaStr1 = dnaStr.upper()
	while i < len(dnaStr1):
		nBytes = nBytes + DNAToByte(dnaStr1[i:i+4])
		i = i + 4
	return nBytes


#single base mapping
def byteToDNA(abyte):
	#define the converter of four bits to DNA chars
	#A=00, C=01, G=10, C=11
	converter = ('A', 'C', 'G', 'T')
	octNum = abyte
	#octNum = ord(abyte)
	dnastr = converter[octNum % 4]
	octNum = octNum//4
	dnastr = converter[octNum % 4] + dnastr
	octNum = octNum // 4
	dnastr = converter[octNum % 4] + dnastr
	octNum = octNum // 4
	dnastr = converter[octNum % 4] + dnastr
	return dnastr



def DNAToByte(dnaStr):
	# converter = ('A', 'T', 'G', 'C')
	#
	# converter = ('A', 'C', 'G', 'T')
	# i = 0
	# while i < 4:
	# 	j = 0
	# 	while j < 4:
	# 		k = 0
	# 		while k < 4:
	# 			m = 0
	# 			while m < 4:
	# 				fourBases = converter[i] + converter[j] + converter[k] + converter[m]
	# 				bytesConverter[fourBases] = bytes([i*4*4*4 + j*4*4 + k*4 + m])
	# 				m = m + 1
	# 			k = k + 1
	# 		j = j + 1
	# 	i = i + 1
	bytesConverter = {'AAAA': b'\x00', 'AAAT': b'\x03', 'AAAG': b'\x02', 'AAAC': b'\x01', 'AATA': b'\x0c', 'AATT': b'\x0f', 'AATG': b'\x0e', 'AATC': b'\r', 'AAGA': b'\x08', 'AAGT': b'\x0b', 'AAGG': b'\n', 'AAGC': b'\t', 'AACA': b'\x04', 'AACT': b'\x07', 'AACG': b'\x06', 'AACC': b'\x05', 'ATAA': b'0', 'ATAT': b'3', 'ATAG': b'2', 'ATAC': b'1', 'ATTA': b'<', 'ATTT': b'?', 'ATTG': b'>', 'ATTC': b'=', 'ATGA': b'8', 'ATGT': b';', 'ATGG': b':', 'ATGC': b'9', 'ATCA': b'4', 'ATCT': b'7', 'ATCG': b'6', 'ATCC': b'5', 'AGAA': b' ', 'AGAT': b'#', 'AGAG': b'"', 'AGAC': b'!', 'AGTA': b',', 'AGTT': b'/', 'AGTG': b'.', 'AGTC': b'-', 'AGGA': b'(', 'AGGT': b'+', 'AGGG': b'*', 'AGGC': b')', 'AGCA': b'$', 'AGCT': b"'", 'AGCG': b'&', 'AGCC': b'%', 'ACAA': b'\x10', 'ACAT': b'\x13', 'ACAG': b'\x12', 'ACAC': b'\x11', 'ACTA': b'\x1c', 'ACTT': b'\x1f', 'ACTG': b'\x1e', 'ACTC': b'\x1d', 'ACGA': b'\x18', 'ACGT': b'\x1b', 'ACGG': b'\x1a', 'ACGC': b'\x19', 'ACCA': b'\x14', 'ACCT': b'\x17', 'ACCG': b'\x16', 'ACCC': b'\x15', 'TAAA': b'\xc0', 'TAAT': b'\xc3', 'TAAG': b'\xc2', 'TAAC': b'\xc1', 'TATA': b'\xcc', 'TATT': b'\xcf', 'TATG': b'\xce', 'TATC': b'\xcd', 'TAGA': b'\xc8', 'TAGT': b'\xcb', 'TAGG': b'\xca', 'TAGC': b'\xc9', 'TACA': b'\xc4', 'TACT': b'\xc7', 'TACG': b'\xc6', 'TACC': b'\xc5', 'TTAA': b'\xf0', 'TTAT': b'\xf3', 'TTAG': b'\xf2', 'TTAC': b'\xf1', 'TTTA': b'\xfc', 'TTTT': b'\xff', 'TTTG': b'\xfe', 'TTTC': b'\xfd', 'TTGA': b'\xf8', 'TTGT': b'\xfb', 'TTGG': b'\xfa', 'TTGC': b'\xf9', 'TTCA': b'\xf4', 'TTCT': b'\xf7', 'TTCG': b'\xf6', 'TTCC': b'\xf5', 'TGAA': b'\xe0', 'TGAT': b'\xe3', 'TGAG': b'\xe2', 'TGAC': b'\xe1', 'TGTA': b'\xec', 'TGTT': b'\xef', 'TGTG': b'\xee', 'TGTC': b'\xed', 'TGGA': b'\xe8', 'TGGT': b'\xeb', 'TGGG': b'\xea', 'TGGC': b'\xe9', 'TGCA': b'\xe4', 'TGCT': b'\xe7', 'TGCG': b'\xe6', 'TGCC': b'\xe5', 'TCAA': b'\xd0', 'TCAT': b'\xd3', 'TCAG': b'\xd2', 'TCAC': b'\xd1', 'TCTA': b'\xdc', 'TCTT': b'\xdf', 'TCTG': b'\xde', 'TCTC': b'\xdd', 'TCGA': b'\xd8', 'TCGT': b'\xdb', 'TCGG': b'\xda', 'TCGC': b'\xd9', 'TCCA': b'\xd4', 'TCCT': b'\xd7', 'TCCG': b'\xd6', 'TCCC': b'\xd5', 'GAAA': b'\x80', 'GAAT': b'\x83', 'GAAG': b'\x82', 'GAAC': b'\x81', 'GATA': b'\x8c', 'GATT': b'\x8f', 'GATG': b'\x8e', 'GATC': b'\x8d', 'GAGA': b'\x88', 'GAGT': b'\x8b', 'GAGG': b'\x8a', 'GAGC': b'\x89', 'GACA': b'\x84', 'GACT': b'\x87', 'GACG': b'\x86', 'GACC': b'\x85', 'GTAA': b'\xb0', 'GTAT': b'\xb3', 'GTAG': b'\xb2', 'GTAC': b'\xb1', 'GTTA': b'\xbc', 'GTTT': b'\xbf', 'GTTG': b'\xbe', 'GTTC': b'\xbd', 'GTGA': b'\xb8', 'GTGT': b'\xbb', 'GTGG': b'\xba', 'GTGC': b'\xb9', 'GTCA': b'\xb4', 'GTCT': b'\xb7', 'GTCG': b'\xb6', 'GTCC': b'\xb5', 'GGAA': b'\xa0', 'GGAT': b'\xa3', 'GGAG': b'\xa2', 'GGAC': b'\xa1', 'GGTA': b'\xac', 'GGTT': b'\xaf', 'GGTG': b'\xae', 'GGTC': b'\xad', 'GGGA': b'\xa8', 'GGGT': b'\xab', 'GGGG': b'\xaa', 'GGGC': b'\xa9', 'GGCA': b'\xa4', 'GGCT': b'\xa7', 'GGCG': b'\xa6', 'GGCC': b'\xa5', 'GCAA': b'\x90', 'GCAT': b'\x93', 'GCAG': b'\x92', 'GCAC': b'\x91', 'GCTA': b'\x9c', 'GCTT': b'\x9f', 'GCTG': b'\x9e', 'GCTC': b'\x9d', 'GCGA': b'\x98', 'GCGT': b'\x9b', 'GCGG': b'\x9a', 'GCGC': b'\x99', 'GCCA': b'\x94', 'GCCT': b'\x97', 'GCCG': b'\x96', 'GCCC': b'\x95', 'CAAA': b'@', 'CAAT': b'C', 'CAAG': b'B', 'CAAC': b'A', 'CATA': b'L', 'CATT': b'O', 'CATG': b'N', 'CATC': b'M', 'CAGA': b'H', 'CAGT': b'K', 'CAGG': b'J', 'CAGC': b'I', 'CACA': b'D', 'CACT': b'G', 'CACG': b'F', 'CACC': b'E', 'CTAA': b'p', 'CTAT': b's', 'CTAG': b'r', 'CTAC': b'q', 'CTTA': b'|', 'CTTT': b'\x7f', 'CTTG': b'~', 'CTTC': b'}', 'CTGA': b'x', 'CTGT': b'{', 'CTGG': b'z', 'CTGC': b'y', 'CTCA': b't', 'CTCT': b'w', 'CTCG': b'v', 'CTCC': b'u', 'CGAA': b'`', 'CGAT': b'c', 'CGAG': b'b', 'CGAC': b'a', 'CGTA': b'l', 'CGTT': b'o', 'CGTG': b'n', 'CGTC': b'm', 'CGGA': b'h', 'CGGT': b'k', 'CGGG': b'j', 'CGGC': b'i', 'CGCA': b'd', 'CGCT': b'g', 'CGCG': b'f', 'CGCC': b'e', 'CCAA': b'P', 'CCAT': b'S', 'CCAG': b'R', 'CCAC': b'Q', 'CCTA': b'\\', 'CCTT': b'_', 'CCTG': b'^', 'CCTC': b']', 'CCGA': b'X', 'CCGT': b'[', 'CCGG': b'Z', 'CCGC': b'Y', 'CCCA': b'T', 'CCCT': b'W', 'CCCG': b'V', 'CCCC': b'U'}
	return bytesConverter[dnaStr.upper()]


#2019-5-14
def calc_gc(dnastr):
	dnastr = dnastr.upper()
	gc_num = dnastr.count("G") + dnastr.count("C")
	return gc_num/len(dnastr)

def max_homo_len(dnastr):
	dnastr = dnastr.upper()
	max_len = 0
	pre_len = 0
	last_c = ''
	for c in dnastr:
		if c == last_c:
			pre_len = pre_len + 1
		else:
			if pre_len > max_len:
				max_len = pre_len
			pre_len = 1
		last_c = c
	if pre_len > max_len:
		max_len = pre_len
	return max_len

#2019-05-15
#2020-03-08
#def check_dna(dnastr, min_gc=0.45, max_gc=0.55, max_homo=5):
def check_dna(dnastr, min_gc=0.4, max_gc=0.6, max_homo=7):
	# print(dnastr)
	gc_rate = calc_gc(dnastr)
	if gc_rate > max_gc:
		return False
	if gc_rate < min_gc:
		return False
	homo_poly_len = max_homo_len(dnastr)
	if homo_poly_len > max_homo:
		return False
	return True

#2019-05-16
def rev_seq(dnastr):
	complement = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}
	dnastr = dnastr.upper()
	rev_dna = ""
	for i in dnastr:
		rev_dna += complement[i]
	rev_dna = rev_dna[::-1]
	return rev_dna

#20190625
def num_randint(a,b,num):
	i = 0
	int_nums = []
	while i < num:
		int_nums.append(random.randint(a, b))
		i += 1
	return int_nums


def randomATGC():
	a = "ATGC"
	ar = random.randint(0, 3)
	return a[ar:ar + 1]


def randomDNA(len):
	i = 1
	dna = ""
	while(i <=len):
		dna= dna + randomATGC()
		i = i + 1

	return dna


def DNA_complement(sequence):
	sequence = sequence.upper()
	sequence = sequence.replace('A', 't')
	sequence = sequence.replace('T', 'a')
	sequence = sequence.replace('C', 'g')
	sequence = sequence.replace('G', 'c')
	return sequence.upper()

def DNA_reverse(sequence):
	sequence = sequence.upper()
	return sequence[::-1]

def DNA_rev_complement(sequence):
	sequence = DNA_complement(sequence)
	return DNA_reverse(sequence)

def expntl(L):
	"""
    negative exponential distribution
    return a double random number, L is the mean value
    """
	u = random.randint(0,2**L)
	return int(np.math.log(u)/math.log(2))


def file_to_array(file):
	f = open(file)
	arr = []
	line = f.readline()
	while line.strip():
		arr.append(line.strip())
		line = f.readline()
	return arr


def kmers_sta_str(str, kmer_len=21,step_len=1):
	kmers = {}
	i = 0
	if len(str) >= kmer_len:
		i = 0
		kmstr = ''
		while i <= len(str) - kmer_len:
			kmstr = str[i:i + kmer_len]
			if kmstr in kmers:
				kmers[kmstr] = kmers[kmstr] + 1
			else:
				kmers[kmstr] = 1
			i = i + step_len
	return kmers

def kmers_of_str(str, kmer_len=21,step_len=1):
	kmers = {}
	i = 0
	if len(str) >= kmer_len:
		i = 0
		kmstr = ''
		while i <= len(str) - kmer_len:
			kmstr = str[i:i + kmer_len]
			kmers[kmstr] = 1
			i = i + step_len
	return kmers.keys()

def kmer_list_of_str(str, kmer_len=21,step_len=1):
	kmers = []
	i = 0
	if len(str) >= kmer_len:
		i = 0
		kmstr = ''
		while i <= len(str) - kmer_len:
			kmstr = str[i:i + kmer_len]
			kmers.append(kmstr)
			i = i + step_len
	return kmers

#2020-05-16
def kmers_of_position(str, kmer_len,pos=0):
	if len(str)-pos >= kmer_len:
		return str[pos:pos+kmer_len]




#2020-05-03
def kmers_in_dict(kmers, dict):
	for kmer in kmers:
		if kmer not in dict:
			return False
	return True

def any_kmers_in_dict(kmers, dict):
	for kmer in kmers:
		if kmer in dict:
			return True
	return False

def read_file(file):
	file1 = open(file, 'rb')
	filebytes = file1.read()
	file1.close()
	return filebytes


def des_en(s, sec_key):
	iv = sec_key
	k = des(sec_key, CBC, iv, pad=None, padmode=PAD_PKCS5)
	return k.encrypt(s)

def des_de(s, sec_key):
	iv = sec_key
	k = des(sec_key, CBC, iv, pad=None, padmode=PAD_PKCS5)
	return k.decrypt(s)

def sha256_rand(seed):
	sd = bytes([seed])
	a = sha256.sha256(sd).digest()
	b = sha256.sha256(a).digest()
	c = sha256.sha256(b).digest()
	return int.from_bytes(c[0:24], 'big')

def sort_bytes_by_order(byte_string, order):
	order_list = [int(x) for x in order]
	new_byte_string = b''
	for i in order_list:
		new_byte_string = new_byte_string + byte_string[i:i + 1]

	return new_byte_string


def reverse_sort_bytes_by_order(byte_string, order):
	order_list = [int(x) for x in order]
	original_byte_string = b''
	for i in range(len(byte_string)):
		original_byte_string += byte_string[order_list.index(i):order_list.index(i) + 1]

	return original_byte_string

def rand_order(max_num, rd_seed):
	list = range(0, max_num)
	random.seed(rd_seed)
	return random.sample(list, max_num)

def ratio(num1, num2):
	if num1 > num2:
		return num1/num2
	else:
		return num2/num1


# def calc_crc8(data):
# 	crc = b'\x00'  # Initial CRC value as a single-byte bytes object
# 	polynomial = b'\x07'  # CRC-8 polynomial (x^8 + x^2 + x^1 + x^0) as a single-byte bytes object
#
# 	for byte in data:
# 		crc = bytes([crc[0] ^ byte])
# 		for _ in range(8):
# 			if crc[0] & 0x80:
# 				crc = bytes([((crc[0] << 1) ^ int.from_bytes(polynomial, 'big'))])
# 			else:
# 				crc = bytes([(crc[0] << 1)])
# 	return crc
#
#
# def is_valid_crc8(data, received_crc):
# 	calculated_crc = calc_crc8(data)
# 	return calculated_crc == received_crc


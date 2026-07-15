import hashlib
import random
import math
import numpy as np
from random import choices

# pyDes is imported lazily by des_en/des_de. Importing it at module scope makes an
# optional, rarely-used cipher a hard requirement for `from utils import *`, which every
# module in this package does -- a missing pyDes would otherwise break the whole codebase.
def _load_des():
	from pyDes import des, CBC, PAD_PKCS5
	return des, CBC, PAD_PKCS5


def	charN(str, N):
	if N < len(str):
		return str[N]
	return 'X'


def oxor(str1, str2):
	length = max(len(str1), len(str2))
	return ''.join(chr(ord(charN(str1, i)) ^ ord(charN(str2, i))) for i in range(length))

def	xor(str1, str2):
	# Zero-pads the shorter operand, then XORs as one big integer instead of
	# concatenating bytes one at a time (which was quadratic).
	length = max(len(str1), len(str2))
	if not length:
		return b''
	a = int.from_bytes(str1, 'big') << (8 * (length - len(str1)))
	b = int.from_bytes(str2, 'big') << (8 * (length - len(str2)))
	return (a ^ b).to_bytes(length, 'big')


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
	# NOTE: this deliberately deviates from Luby's ideal soliton distribution, which
	# defines rho(1) = 1/K. Using 10/K puts ~5x more mass on degree 1 (after
	# normalisation) and range(2, K) means degree K is never emitted. The result is
	# still a valid normalised distribution. Preserved as-is so that previously tuned
	# experimental results stay reproducible -- see README "Known deviations".
	if K < 2:
		raise ValueError("ideal_distribution needs K >= 2, got %r" % (K,))

	probabilities = [0, 10 / K]
	probabilities += [1 / (k * (k - 1)) for k in range(2, K)]
	probabilities_sum = np.sum(probabilities)
	probabilities /= probabilities_sum

	# assert probabilities_sum >= 1 - epsilon and probabilities_sum <= 1 + epsilon, "The ideal distribution should be standardized"
	return probabilities

def robust_distribution(K, c_value = 0.5, robust_failure_probability=0.1):
	# NOTE: the tau spike lands at index M-1 (degree M-1) rather than Luby's degree M,
	# and tau covers degrees 1..M-2 rather than 1..M-1. Preserved as-is for
	# reproducibility -- see README "Known deviations".
	if K < 2:
		raise ValueError("robust_distribution needs K >= 2, got %r" % (K,))

	S = c_value * math.log(K / robust_failure_probability) * math.sqrt(K)
	M = round(K / S)

	# M > K used to build an over-long extra_proba and die in np.add with an opaque
	# broadcast error; M < 2 produced an empty tau. Clamp so the spike stays in range.
	M = min(max(M, 2), K)

	extra_proba = [0] + [1/(i * M) for i in range(1, M-1)]
	extra_proba += [S * math.log(S / robust_failure_probability) / K]  # Spike at M
	extra_proba += [0 for k in range(M, K)]

	probabilities = np.add(extra_proba, ideal_distribution(K))
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
	return ''.join(map(_BYTE_TO_DNA.__getitem__, manyBytes))


def DNAToBytes(dnaStr):
	# Vectorised 2-bit packing. Equivalent to repeatedly calling DNAToByte on 4-base
	# groups, but ~2400x faster on genome-scale input (12.7s -> 5ms for E. coli 4.6Mbp).
	if not dnaStr:
		return b''
	try:
		codes = dnaStr.encode('ascii').translate(_DNA_TO_2BIT)
	except UnicodeEncodeError:
		raise ValueError("DNA sequence contains non-ASCII characters")
	if len(codes) % 4:
		raise ValueError(
			"DNA length must be a multiple of 4, got %d (pad with ChromosomeECC.round_DNA)"
			% len(codes))
	arr = np.frombuffer(codes, dtype=np.uint8)
	if arr.max() > 3:
		bad = dnaStr[int(np.argmax(arr > 3))]
		raise ValueError("DNA sequence contains non-ACGT base %r" % bad)
	arr = arr.reshape(-1, 4)
	packed = (arr[:, 0] << 6) | (arr[:, 1] << 4) | (arr[:, 2] << 2) | arr[:, 3]
	return packed.astype(np.uint8).tobytes()


#single base mapping
#A=00, C=01, G=10, T=11
_BASES = ('A', 'C', 'G', 'T')

# byte -> its 4-base string, and the reverse. Built once at import; the previous code
# rebuilt a 256-entry dict on every single 4-base conversion.
_BYTE_TO_DNA = tuple(
	_BASES[(b >> 6) & 3] + _BASES[(b >> 4) & 3] + _BASES[(b >> 2) & 3] + _BASES[b & 3]
	for b in range(256))
_DNA_TO_BYTE = {seq: bytes([b]) for b, seq in enumerate(_BYTE_TO_DNA)}

# ASCII base -> 2-bit code; every non-ACGT byte maps to 255 so it fails the range check.
_DNA_TO_2BIT = bytes(
	'ACGT'.index(chr(i).upper()) if chr(i).upper() in 'ACGT' else 255
	for i in range(256))


def byteToDNA(abyte):
	return _BYTE_TO_DNA[abyte]



def DNAToByte(dnaStr):
	# Raises KeyError for non-ACGT input or a group that is not exactly 4 bases,
	# matching the original hand-written lookup table this replaced.
	return _DNA_TO_BYTE[dnaStr.upper()]


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
	# np.math was removed in numpy 2.0, and randint's lower bound of 0 made log(0)
	# a domain error roughly once every 2**L draws.
	u = random.randint(1, 2**L)
	return int(math.log(u)/math.log(2))


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


# WARNING: DES has a 56-bit key and is not considered secure. These helpers are kept
# only for compatibility with data encoded by earlier versions of this codebase, and
# the IV is the key itself, which makes CBC deterministic. Do not use for secrecy.
# See README "Security notes".
def des_en(s, sec_key):
	des, CBC, PAD_PKCS5 = _load_des()
	iv = sec_key
	k = des(sec_key, CBC, iv, pad=None, padmode=PAD_PKCS5)
	return k.encrypt(s)

def des_de(s, sec_key):
	des, CBC, PAD_PKCS5 = _load_des()
	iv = sec_key
	k = des(sec_key, CBC, iv, pad=None, padmode=PAD_PKCS5)
	return k.decrypt(s)

def sha256_rand(seed):
	# Was bytes([seed]), which raised ValueError for any seed outside 0..255 and gave
	# only 256 distinct streams. Uses stdlib hashlib rather than the third-party
	# 'sha256' package, which reimplemented the same primitive.
	sd = int(seed).to_bytes(max(1, (int(seed).bit_length() + 7) // 8), 'big', signed=False)
	a = hashlib.sha256(sd).digest()
	b = hashlib.sha256(a).digest()
	c = hashlib.sha256(b).digest()
	return int.from_bytes(c[0:24], 'big')

def sort_bytes_by_order(byte_string, order):
	order_list = [int(x) for x in order]
	return bytes(byte_string[i] for i in order_list)


def reverse_sort_bytes_by_order(byte_string, order):
	# Inverts sort_bytes_by_order. The previous version called order_list.index(i)
	# inside the loop, making it quadratic; a single inverse-permutation pass is O(n).
	order_list = [int(x) for x in order]
	original = bytearray(len(byte_string))
	for src, dst in enumerate(order_list):
		if src < len(byte_string):
			original[dst] = byte_string[src]
	return bytes(original)

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


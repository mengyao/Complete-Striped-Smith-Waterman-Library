from ctypes import *


# ctypes Structure with s_align struct mapping
# This structure is returned by SSWAligner.Align func
class SSWAlignRes(Structure):
		_fields_ = [('score', c_uint16),
			('score2', c_uint16),
			('ref_begin1', c_int32),
			('ref_end1', c_int32),
			('read_begin2', c_int32),
			('read_end2', c_int32),
			('ref_end2', c_int32),
			('cigar', POINTER(c_uint32)),
			('cigarLen', c_int32)]

# Wrapper for SSW align library
class SSWAligner:

	# Hashmap to map Nucleotide to int expected by the SSW C library
	base_to_int = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 , 'N': 4}
	
	
	def __init__(self, read_seq, match = 2, mismatch = -2):
		# Load the ssw library using ctypes
		libssw = cdll.LoadLibrary('libssw.so')
		
		# Init and setup the functions pointer to latch the one specified in the SSW lib 
		self.ssw_init = libssw.ssw_init
		self.ssw_init.restype = c_void_p
		self.ssw_init.argtypes = [POINTER(c_int8), c_int32, POINTER(c_int8), c_int32, c_int8]

		self.init_destroy = libssw.init_destroy
		self.init_destroy.restype = None
		self.init_destroy.argtypes =  [c_void_p]

		self.ssw_align = libssw.ssw_align
		self.ssw_align.restype = POINTER(SSWAlignRes)
		self.ssw_align.argtypes = [c_void_p, POINTER(c_int8), c_int32, c_uint8, c_uint8, c_uint8, c_uint16, c_int32, c_int32]

		self.align_destroy = libssw.align_destroy
		self.align_destroy.restype = None
		self.align_destroy.argtypes = [POINTER(SSWAlignRes)]

		self.read_seq = read_seq

		self.read_len = len(self.read_seq)	
		# Cost matrix 	
		mat_decl = c_int8 * 25
		self.mat = mat_decl(match, mismatch, mismatch, mismatch, 0,
				mismatch, match, mismatch, mismatch, 0,
				mismatch, mismatch, match, mismatch, 0,
				mismatch, mismatch, mismatch, match, 0,
				 0, 0, 0, 0, 0)
		
		read_num_decl = c_int8 * self.read_len
		self.read_num = read_num_decl()
		k = 0
		for i in self.read_seq:
			value = self.base_to_int.get(i)
			if value == None:
				value = 4
			self.read_num[k] = value 
			k += 1
		
		self.profile = self.ssw_init(self.read_num, c_int32(self.read_len), self.mat, 5, 2) 
	
	def Destroy(self):
		self.init_destroy(self.profile)

	def Align(self, str_ref, gap_open = 3, gap_extension = 1):
		str_ref_len = len(str_ref)
		ref_num_decl= c_int8 * str_ref_len
		ref_num = ref_num_decl()
		k = 0
		for i in str_ref:
			value = self.base_to_int.get(i)
			if value == None:
				value = 4
				
			ref_num[k] = value
			k += 1
		
		mask_len = self.read_len / 2
		if mask_len < 15:
			mask_len = 15
	
		return self.ssw_align(self.profile, ref_num, c_int32(str_ref_len), gap_open, gap_extension, 1, 0, 0, mask_len) 


	def DestroyAlignRes(self, align):
		self.align_destroy(align)


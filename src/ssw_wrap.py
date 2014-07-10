"""
@package    ssw_wrap
@brief      Simple python wrapper for SSW align library
@copyright  [The MIT licence](http://opensource.org/licenses/MIT)
@author     Clement & Adrien Leger - 2014
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages
from ctypes import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Aligner(object):
    """
    @class  SSWAligner
    @brief Wrapper for SSW align library
    To use the dynamic library libssw.so you may need to modify the LD_LIBRARY_PATH environment
    variable to include the library directory
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS VARIABLES~~~~~~~#

    # Dictionnary to map Nucleotide to int as expected by the SSW C library
    base_to_int = { 'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
    int_to_base = { 0:'A', 1:'C', 2:'G', 3:'T', 4:'N'}

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "Match bonus\t{}\n".format(self.match)
        msg += "Mismatch penalty\t{}\n".format(self.mismatch)
        msg += "Gap open penalty\t{}\n".format(self.gap_open)
        msg += "Gap extend penalty\t{}\n".format(self.gap_open)
        msg += "Score matrix\n"
        msg += "\tA\tC\tG\tT\tN\n"
        msg += "A\t{}\t{}\t{}\t{}\t{}\n".format(self.match, -self.mismatch, -self.mismatch, -self.mismatch, 0)
        msg += "C\t{}\t{}\t{}\t{}\t{}\n".format(-self.mismatch, self.match, -self.mismatch, -self.mismatch, 0)
        msg += "G\t{}\t{}\t{}\t{}\t{}\n".format(-self.mismatch, -self.mismatch, self.match, -self.mismatch, 0)
        msg += "T\t{}\t{}\t{}\t{}\t{}\n".format(-self.mismatch, -self.mismatch, -self.mismatch, self.match, 0)
        msg += "N\t{}\t{}\t{}\t{}\t{}\n".format(0,0,0,0,0)
        msg += "Reference Sequence :\n"
        msg += "".join([self.int_to_base[i] for i in self.ref_seq])+"\n"
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__(self, ref_seq, match=2, mismatch=2, gap_open=3, gap_extension=1):
        """
        Initialize object by creating an interface with ssw library fonctions
        A reference sequence is also assigned to the object for multiple alignment against queries
        with the align function
        @param ref_seq Reference sequence as a python string
        @param match Weight for a match
        @param mismatch Absolute value of mismatch penalty
        @param gap_open Absolute value of gap open penalty
        @param gap_extension Absolute value of gap extend penalty
        """
        # Store object variables
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extension = gap_extension

        # Load the ssw library using ctypes
        libssw = cdll.LoadLibrary('libssw.so')

        # Init and setup the functions pointer to latch the one specified in the SSW lib
        # ssw_init method
        self.ssw_init = libssw.ssw_init
        self.ssw_init.restype = c_void_p
        self.ssw_init.argtypes = [POINTER(c_int8), c_int32, POINTER(c_int8), c_int32, c_int8]

        # init_destroy method
        self.init_destroy = libssw.init_destroy
        self.init_destroy.restype = None
        self.init_destroy.argtypes =  [c_void_p]

        # ssw_align method
        self.ssw_align = libssw.ssw_align
        self.ssw_align.restype = POINTER(CAlignRes)
        self.ssw_align.argtypes = [c_void_p, POINTER(c_int8), c_int32, c_uint8, c_uint8, c_uint8, c_uint16, c_int32, c_int32]

        # align_destroy method
        self.align_destroy = libssw.align_destroy
        self.align_destroy.restype = None
        self.align_destroy.argtypes = [POINTER(CAlignRes)]

        # Initialize Cost matrix . ambiguous base: no penalty
        mat_decl = c_int8 * 25
        self.mat = mat_decl(match, -mismatch, -mismatch, -mismatch, 0,
                            -mismatch, match, -mismatch, -mismatch, 0,
                            -mismatch, -mismatch, match, -mismatch, 0,
                            -mismatch, -mismatch, -mismatch, match, 0,
                            0, 0, 0, 0, 0)

        # Determine the size of the ref sequence and cast it in a c type integer matrix
        self.ref_len = len(ref_seq)
        self.ref_seq = self._DNA_to_int_mat (ref_seq, self.ref_len)


    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def align(self, query_seq):
        """
        Perform the alignment of query against the object reference sequence
        @param query_seq Query sequence as a python string
        @return A SSWAlignRes Object containing informations about the alignment.
        """
        # Determine the size of the ref sequence and cast it in a c type integer matrix
        query_len = len(query_seq)
        query_seq = self._DNA_to_int_mat (query_seq, query_len)

        # Create the query profile using the query sequence
        profile = self.ssw_init(query_seq, # Query seq in c type integers
                                c_int32(query_len), # Length of Queryseq in bites
                                self.mat, # Score matrix
                                5, # Square root of the number of elements in mat
                                2) # flag = no estimation of the best alignment score

        # Setup the mask_len parameters = distance between the optimal and suboptimal alignment
        # if < 15, the function will NOT return the suboptimal alignment information
        mask_len = query_len / 2
        if mask_len < 15:
            mask_len = 15

        c_result = self.ssw_align (profile, # Query profile
                                self.ref_seq, # Ref seq in c type integers
                                c_int32(self.ref_len), # Length of Refseq in bites
                                self.gap_open, # Absolute value of gap open penalty
                                self.gap_extension, # absolute value of gap extend penalty
                                1, # Bitwise FLAG for output values = not set
                                0, # Score filter
                                0, # Distance filter
                                mask_len) # Distance between the optimal and suboptimal alignment

        return c_result
        ## Transform the Cstructure into a python object
        #py_result = PyAlignRes(c_result)

        ## Free reserved space by ssw.init and ssw_init methods.
        #self._init_destroy(profile)
        #self._align_destroy(c_result)

        ## Return the object
        #return py_result

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _DNA_to_int_mat (self, seq, len_seq):
        """
        Cast a python DNA string into a Ctype int8 matrix
        """

        # Declare the matrix
        read_num_decl = c_int8 * len_seq
        read_num = read_num_decl()

        # for each letters in ATCGN transform in integers thanks to self.base_to_int
        for i in range(len_seq):
            try:
                value = self.base_to_int[seq[i]]
            # if the base is not in the canonic DNA bases assign 4 as for N
            except KeyError:
                value = 4
            finally:
                read_num[i] = value

        return read_num

    def _init_destroy(self, profile):
        """
        Free the space alocated for the matrix used by init
        """
        self.init_destroy(profile)

    def _align_destroy(self, align):
        """
        Free the space alocated for the matrix used by align
        """
        self.align_destroy(align)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class CAlignRes(Structure):
    """
    @class  SSWAlignRes
    @brief  ctypes Structure with s_align struct mapping returned by SSWAligner.Align func
            Correspond to the structure of the query profile
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~Ctype Structure~~~~~~~#
    _fields_ = [('score', c_uint16),
                ('score2', c_uint16),
                ('ref_begin1', c_int32),
                ('ref_end1', c_int32),
                ('read_begin2', c_int32),
                ('read_end2', c_int32),
                ('ref_end2', c_int32),
                ('cigar', POINTER(c_uint32)),
                ('cigarLen', c_int32)]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class PyAlignRes(object):
    """
    @class  SSWAlignRes
    @brief  ctypes Structure with s_align struct mapping returned by SSWAligner.Align func
            Correspond to the structure of the query profile
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHOD~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        if self.score1:
            msg += "Score 1\t{}\n".format(self.score1)
        if self.score2:
            msg += "Score 2\t{}\n".format(self.score2)
        if self.cigar:
            msg += "Cigar\t{}\n".format(self.cigar)
        if self.cigarLen:
            msg += "CigarLen\t{}\n".format(self.cigarLen)
        if self.ref_begin1:
            msg += "Ref_begin1\t{}\n".format(self.ref_begin1)
        if self.ref_end1:
            msg += "Ref_end1\t{}\n".format(self.ref_end1)
        if self.ref_begin2:
            msg += "Ref_begin2\t{}\n".format(self.ref_begin2)
        if self.ref_end2:
            msg += "Ref_end2\t{}\n".format(self.ref_end2)
        if self.read_begin1:
            msg += "Read_begin1\t{}\n".format(self.read_begin1)
        if self.read_end1:
            msg += "Read_end1\t{}\n".format(self.read_end1)
        if self.read_begin2:
            msg += "Read_begin2\t{}\n".format(self.read_begin2)
        if self.read_end2:
            msg += "Read_end2\t{}\n".format(self.read_end2)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__ (self, Res):
        """
        """
         # Parse value in the C type structure pointer
        try:
            self.score1 = Res.contents.score
        except AttributeError:
            self.score1 = None
        try:
            self.score2 = Res.contents.score2
        except AttributeError:
            self.score2 = None

        #TODO = write a cigar as a string
        try:
            self.cigar = Res.contents.cigar.contents.value
        except AttributeError:
            self.cigar = None
        try:
            self.cigarLen = Res.contents.cigarLen
        except AttributeError:
            self.cigarLen = None
        try:
            self.ref_begin1 = Res.contents.ref_begin1
        except AttributeError:
            self.ref_begin1 = None
        try:
            self.ref_begin2 = Res.contents.ref_begin2
        except AttributeError:
            self.ref_begin2 = None
        try:
            self.ref_end1 = Res.contents.ref_end1
        except AttributeError:
            self.ref_end1 = None
        try:
            self.ref_end2 = Res.contents.ref_end2
        except AttributeError:
            self.ref_end2 = None
        try:
            self.read_begin1 = Res.contents.read_begin1
        except AttributeError:
            self.read_begin1 = None
        try:
            self.read_begin2 = Res.contents.read_begin2
        except AttributeError:
            self.read_begin2 = None
        try:
            self.read_end1 = Res.contents.read_end1
        except AttributeError:
            self.read_end1 = None
        try:
            self.read_end2 = Res.contents.read_end2
        except AttributeError:
            self.read_end2 = None

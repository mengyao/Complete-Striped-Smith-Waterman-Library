"""
@package    SswWrapper
@brief      Simple python wrapper for SSW align library
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Clement & Adrien Leger - 2014
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages
from ctypes import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class AlignRes(Structure):
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
class Aligner(object):
    """
    @class  SSWAligner
    @brief Wrapper for SSW align library
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS VARIABLES~~~~~~~#

    # Dictionnary to map Nucleotide to int as expected by the SSW C library
    base_to_int = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 , 'N': 4}

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()

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
        self.ssw_align.restype = POINTER(SSWAlignRes)
        self.ssw_align.argtypes = [c_void_p, POINTER(c_int8), c_int32, c_uint8, c_uint8, c_uint8, c_uint16, c_int32, c_int32]

        # align_destroy method
        self.align_destroy = libssw.align_destroy
        self.align_destroy.restype = None
        self.align_destroy.argtypes = [POINTER(SSWAlignRes)]

        # Initialize Cost matrix . ambiguous base: no penalty
        mat_decl = c_int8 * 25
        self.mat = mat_decl(match, -mismatch, -mismatch, -mismatch, 0,
                            -mismatch, match, -mismatch, -mismatch, 0,
                            -mismatch, -mismatch, match, -mismatch, 0,
                            -mismatch, -mismatch, -mismatch, match, 0,
                            0, 0, 0, 0, 0)

        # Determine the size of the ref sequence and cast it in a c type integer matrix
        self.ref_len = len(ref_seq)
        self.ref_seq = self._DNA_to_int_mat (ref_seq, ref_len)


    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def align(self, query_seq):
        """
        Perform the alignment of query against the object reference sequence
        @param query_seq Query sequence as a python string
        @return A SSWAlignRes Object containing informations about the alignment.
        """
        # Determine the size of the ref sequence and cast it in a c type integer matrix
        query_len = len(query_seq)
        query_seq = self._DNA_to_int_mat (query_seq, ref_len)

        # Create the query profile using the query sequence
        q_profile = self.ssw_init(query_seq, # Query seq in c type integers
                                c_int32(query_len), # Length of Queryseq in bites
                                self.mat, # Score matrix
                                5, # Square root of the number of elements in mat
                                2) # flag = no estimation of the best alignment score

        # Setup the mask_len parameters = distance between the optimal and suboptimal alignment
        # if < 15, the function will NOT return the suboptimal alignment information
        mask_len = query_len / 2
        if mask_len < 15:
            mask_len = 15

        result = self.ssw_align (profile, # Query profile
                                self.ref_seq, # Ref seq in c type integers
                                c_int32(self.ref_len), # Length of Refseq in bites
                                self.gap_open, # Absolute value of gap open penalty
                                self.gap_extension, # absolute value of gap extend penalty
                                0, # Bitwise FLAG for output values = not set
                                0, # Score filter
                                0, # Distance filter
                                mask_len) # Distance between the optimal and suboptimal alignment

        # Free reserved space by ssw.init and ssw_init methods.
        self._init_destroy()
        self._align_destroy()

        return result

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

    def _init_destroy(self):
        """
        Free the space alocated for the matrix used by init
        """
        self.init_destroy(self.profile)

    def _align_destroy(self, align):
        """
        Free the space alocated for the matrix used by align
        """
        self.align_destroy(align)

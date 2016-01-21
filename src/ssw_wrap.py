"""
@package ssw_wrap
@brief Simple python wrapper for SSW align library
To use the dynamic library libssw.so you may need to modify the LD_LIBRARY_PATH environment
variable to include the library directory (export LD_LIBRARY_PATH=$PWD) or for definitive
inclusion of the lib edit /etc/ld.so.conf and add the path or the directory containing the
library and update the cache by using /sbin/ldconfig as root
@copyright  [The MIT licence](http://opensource.org/licenses/MIT)
@author     Clement & Adrien Leger - 2014
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages
import os
import itertools
from ctypes import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

def load_ssw_library():
    relative_path = os.path.split(__file__)[0]
    paths = [".", ".."]
    libname = "libssw.so"
    for path in paths:
        libpath = os.path.join(relative_path, path, libname)
        try:
            return cdll.LoadLibrary(libpath)
        except OSError:
            pass
    # last attempt
    return cdll.LoadLibrary(libname)

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
                ('ref_begin', c_int32),
                ('ref_end', c_int32),
                ('query_begin', c_int32),
                ('query_end', c_int32),
                ('ref_end2', c_int32),
                ('cigar', POINTER(c_uint32)),
                ('cigarLen', c_int32)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Aligner(object):
    """
    @class  SSWAligner
    @brief Wrapper for SSW align library
    """

    # Dictionnary to map Nucleotide to int as expected by the SSW C library
    base_to_int = { 'A':0, 'C':1, 'G':2, 'T':3, 'N':4, 'a':0, 'c':1, 'g':2, 't':3, 'n':4}
    int_to_base = { 0:'A', 1:'C', 2:'G', 3:'T', 4:'N'}

    # Load the ssw library using ctypes
    libssw = load_ssw_library()

    # Init and setup the functions pointer to map the one specified in the SSW lib
    # ssw_init method
    ssw_init = libssw.ssw_init
    ssw_init.restype = c_void_p
    ssw_init.argtypes = [POINTER(c_int8), c_int32, POINTER(c_int8), c_int32, c_int8]
    # init_destroy function
    init_destroy = libssw.init_destroy
    init_destroy.restype = None
    init_destroy.argtypes =  [c_void_p]
    # ssw_align function
    ssw_align = libssw.ssw_align
    ssw_align.restype = POINTER(CAlignRes)
    ssw_align.argtypes = [c_void_p, POINTER(c_int8), c_int32, c_uint8, c_uint8, c_uint8, c_uint16, c_int32, c_int32]
    # align_destroy function
    align_destroy = libssw.align_destroy
    align_destroy.restype = None
    align_destroy.argtypes = [POINTER(CAlignRes)]

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "SCORE PARAMETERS:\n"
        msg += " Gap Weight     Open: {}     Extension: {}\n".format(-self.gap_open, -self.gap_extend)
        msg += " Align Weight   Match: {}    Mismatch: {}\n\n".format(self.match, -self.mismatch)
        msg += " Match/mismatch Score matrix\n"
        msg += " \tA\tC\tG\tT\tN\n"
        msg += " A\t{}\t{}\t{}\t{}\t{}\n".format(self.match, -self.mismatch, -self.mismatch, -self.mismatch, 0)
        msg += " C\t{}\t{}\t{}\t{}\t{}\n".format(-self.mismatch, self.match, -self.mismatch, -self.mismatch, 0)
        msg += " G\t{}\t{}\t{}\t{}\t{}\n".format(-self.mismatch, -self.mismatch, self.match, -self.mismatch, 0)
        msg += " T\t{}\t{}\t{}\t{}\t{}\n".format(-self.mismatch, -self.mismatch, -self.mismatch, self.match, 0)
        msg += " N\t{}\t{}\t{}\t{}\t{}\n\n".format(0,0,0,0,0)
        msg += "REFERENCE SEQUENCE :\n"
        if self.ref_len <= 50:
            msg += "".join([self.int_to_base[i] for i in self.ref_seq])+"\n"
        else:
            msg += "".join([self.int_to_base[self.ref_seq[i]] for i in range(50)])+"...\n"
        msg += " Lenght :{} nucleotides\n".format(self.ref_len)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__(self,
                ref_seq="",
                match=2,
                mismatch=2,
                gap_open=3,
                gap_extend=1,
                report_secondary=False,
                report_cigar=False):
        """
        Initialize object by creating an interface with ssw library fonctions
        A reference sequence is also assigned to the object for multiple alignment against queries
        with the align function
        @param ref_seq Reference sequence as a python string (case insensitive)
        @param match Weight for a match
        @param mismatch Absolute value of mismatch penalty
        @param gap_open Absolute value of gap open penalty
        @param gap_extend Absolute value of gap extend penalty
        @param report_secondary Report the 2nd best alignement if true
        @param report_cigar Report cigar string if true
        """

        # Store overall alignment parameters
        self.report_secondary = report_secondary
        self.report_cigar = report_cigar

        # Set gap penalties
        self.set_gap(gap_open, gap_extend)

        # Set the cost matrix
        self.set_mat(match, mismatch)

        # Set the reference sequence
        self.reference = ref_seq

    #~~~~~~~SETTERS METHODS~~~~~~~#

    def set_gap(self, gap_open=3, gap_extend=1):
        """
        Store gapopen and gap extension penalties
        """
        self.gap_open = gap_open
        self.gap_extend = gap_extend


    def set_mat(self, match=2, mismatch=2):
        """
        Store match and mismatch scores then initialize a Cost matrix and fill it with match and
        mismatch values. Ambiguous base: no penalty
        """
        self.match = match
        self.mismatch = mismatch

        mat_decl = c_int8 * 25
        self.mat = mat_decl(match, -mismatch, -mismatch, -mismatch, 0,
                            -mismatch, match, -mismatch, -mismatch, 0,
                            -mismatch, -mismatch, match, -mismatch, 0,
                            -mismatch, -mismatch, -mismatch, match, 0,
                            0, 0, 0, 0, 0)

    def get_reference(self, ref_seq):
        return self.ref_seq

    def set_reference(self, ref_seq):
        self.ref_seq = ref_seq
        self._ref_seq = self._DNA_to_int_mat(self.ref_seq)
    reference = property(get_reference, set_reference)

    def align(self, query_seq, min_score=0, min_len=0):
        """
        Perform the alignment of query against the object reference sequence
        @param query_seq Query sequence as a python string (case insensitive)
        @param min_score Minimal score of match. None will be return in case of filtering out
        @param min_len Minimal length of match. None will be return in case of filtering out
        @return A SSWAlignRes Object containing informations about the alignment.
        """
        # Determine the size of the ref sequence and cast it in a c type integer matrix
        _query_seq = self._DNA_to_int_mat(query_seq)

        # Create the query profile using the query sequence
        profile = self.ssw_init(_query_seq, # Query seq in c type integers
                                c_int32(len(query_seq)), # Length of Queryseq in bites
                                self.mat, # Score matrix
                                5, # Square root of the number of elements in mat
                                2) # flag = no estimation of the best alignment score

        # Setup the mask_len parameters = distance between the optimal and suboptimal alignment
        # if < 15, the function will NOT return the suboptimal alignment information

        if len(query_seq) > 30:
            mask_len = len(query_seq) / 2
        else:
            mask_len = 15

        c_result = self.ssw_align (profile, # Query profile
                                self._ref_seq, # Ref seq in c type integers
                                c_int32(len(self.ref_seq)), # Length of Refseq in bites
                                self.gap_open, # Absolute value of gap open penalty
                                self.gap_extend, # absolute value of gap extend penalty
                                1, # Bitwise FLAG for output values = return all
                                0, # Score filter = return all
                                0, # Distance filter = return all
                                mask_len) # Distance between the optimal and suboptimal alignment

        # Transform the Cstructure into a python object if score and lenght match the requirements
        score = c_result.contents.score
        match_len  = c_result.contents.query_end - c_result.contents.query_begin + 1

        if score >= min_score and match_len >= min_len:
            py_result = PyAlignRes(c_result, query_seq, self.ref_seq)
        else:
            py_result = None

        # Free reserved space by ssw.init and ssw_init methods.
        self._init_destroy(profile)
        self._align_destroy(c_result)

        # Return the object
        return py_result

    def _DNA_to_int_mat(self, seq):
        # Declare the matrix
        query_num_decl = c_int8 * len(seq)
        query_num = query_num_decl()

        # for each letters in ATCGN transform in integers thanks to self.base_to_int
        for (idx, base) in enumerate(seq):
            try:
                value = self.base_to_int[base]
            # if the base is not in the canonic DNA bases assign 4 as for N
            except KeyError:
                value = 4
            finally:
                query_num[idx] = value

        return query_num

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
class PyAlignRes(object):
    """
    @class  PyAlignRes
    @brief  Extract and verify result from a CAlignRes structure. A comprehensive python
    object is created according to user requirements (+- cigar string and secondary alignment)
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS VARIABLES~~~~~~~#

    # Load the ssw library using ctypes
    libssw = load_ssw_library()

    # Init and setup the functions pointer to map the one specified in the SSW lib
    # cigar_int_to_len function
    cigar_int_to_len = libssw.cigar_int_to_len
    cigar_int_to_len.restype = c_int32
    cigar_int_to_len.argtypes = [c_int32]
    # cigar_int_to_op function
    cigar_int_to_op = libssw.cigar_int_to_op
    cigar_int_to_op.restype = c_char
    cigar_int_to_op.argtypes = [c_int32]

    #~~~~~~~FONDAMENTAL METHOD~~~~~~~#

    def __str__(self):
        msg += "OPTIMAL MATCH\n"
        msg += "Score            {}\n".format(self.score)
        msg += "Reference begin  {}\n".format(self.ref_begin)
        msg += "Reference end    {}\n".format(self.ref_end)
        msg += "Query begin      {}\n".format(self.query_begin)
        msg += "Query end        {}\n".format(self.query_end)

        if self.cigar_string:
            msg += "Cigar_string     {}\n".format(self.cigar_string)

        if self.score2:
            msg += "SUB-OPTIMAL MATCH\n"
            msg += "Score 2           {}\n".format(self.score2)
            msg += "Ref_end2          {}\n".format(self.ref_end2)

        return msg

    def __init__ (self, Res, query_seq, ref_seq):
        self.score = Res.contents.score
        self.ref_seq = ref_seq
        self.ref_begin = Res.contents.ref_begin
        self.ref_end = Res.contents.ref_end
        self.query_seq = query_seq
        self.query_begin = Res.contents.query_begin
        self.query_end = Res.contents.query_end
        self._cigar_string = [Res.contents.cigar[idx] for idx in range(Res.contents.cigarLen)]

    @property
    def iter_cigar(self):
        for val in self._cigar_string:
            op_len = self.cigar_int_to_len(val)
            op_char = self.cigar_int_to_op(val)
            yield (op_len, op_char)

    @property
    def cigar_string(self):
        # Empty string for iterative writing of the cigar string
        cigar_string = ""
        if len(self._cigar_string) == 0:
            return cigar_string

        # If the query match do not start at its first base
        # = introduce a softclip at the begining
        if self.query_begin > 0:
            op_len = self.query_begin
            op_char = "S"
            cigar_string += '{}{}'.format(op_len, op_char)

        # Iterate over the cigar (pointer to a vector of int)
        cigar_string += str.join('', [str.join('', map(str, cstr)) for cstr in self.iter_cigar])

        # If the lenght of bases aligned is shorter than the overall query length
        # = introduce a softclip at the end
        end_len = len(self.query_seq) - self.query_end - 1
        if  end_len != 0:
            op_len = end_len
            op_char = "S"
            cigar_string += '{}{}'.format(op_len, op_char)

        return cigar_string
    cigar = cigar_string

    @property
    def alignment(self):
        def seqiter(seq):
            seq = iter(seq)
            def getseq(cnt):
                return str.join('', [seq.next() for x in range(cnt)])
            return getseq
        print self.ref_seq, self.query_seq
        r_seq = seqiter(self.ref_seq)
        q_seq = seqiter(self.query_seq)
        r_line = m_line = q_line = ''
        for (op_len, op_char) in self.iter_cigar:
            op_len = int(op_len)
            if op_char.upper() == 'M':
                for (r_base, q_base) in zip(r_seq(op_len), q_seq(op_len)):
                    r_line += r_base
                    q_line += q_base
                    if r_base == q_base:
                        m_line += '|'
                    else:
                        m_line += '*'
            elif op_char.upper() == 'I':
                r_line += ' ' * op_len
                m_line += ' ' * op_len
                q_line += q_seq(op_len)
            elif op_char.upper() == 'D':
                r_line += r_seq(op_len)
                m_line += ' ' * op_len
                q_line += ' ' * op_len
        return (r_line, m_line, q_line)

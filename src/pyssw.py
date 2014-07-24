#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@package pyssw
@brief Python standalone program for ssw alignment using the C library
Complete-Striped-Smith-Waterman-Library
Biopython module is require for fastq/fastq parsing
@copyright  [The MIT licence](http://opensource.org/licenses/MIT)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages
import optparse
from imp import find_module
import sys
from time import time

#~~~~~~~MAIN FUNCTION~~~~~~~#
def main (opt):
    
    start = time()
    # Import fasta subject
    if opt.subject.rpartition(".")[2].lower() in ["fa", "fasta"]:
        subject = SeqIO.read(opt.subject, "fasta")
    else:
        raise Exception ("The subject sequence if not in fasta format")
    
    # Import fasta or fastq query
    if opt.query.rpartition(".")[2].lower() in ["fa", "fasta"]:
        query_gen = SeqIO.parse(opt.query, "fasta")
    elif opt.query.rpartition(".")[2].lower() == "fastq":
        query_gen = SeqIO.parse(opt.query, "fastq")
    else:
        raise Exception ("The query sequence if not in fastq or fasta format")
    
    print ("Initialize ssw aligner with the subject sequence")
    # Init the an Aligner object with the reference value 
    ssw = Aligner(
        str(subject.seq),
        match=int(opt.match),
        mismatch=int(opt.mismatch),
        gap_open=int(opt.gap_open),
        gap_extend= int(opt.gap_extend),
        report_secondary=False,
        report_cigar=True)
    
    # Write the header of the SAM file
    with open("result.sam", "w") as f:
        f.write("@HD\tVN:1.0\tSO:unsorted\n")
        f.write("@SQ\tSN:{}\tLN:{}\n".format(subject.id, len(subject.seq)))
        f.write("@PG\tID:Striped-Smith-Waterman\tPN:pyssw\tVN:0.1\n")
        f.write("@CO\tScore_values = match {}, mismatch {}, gap_open {}, gap_extend {}\n".format(
            opt.match,
            opt.mismatch,
            opt.gap_open,
            opt.gap_extend))
        f.write("@CO\tFilter Options = min_score {}, min_len {}\n".format(
            opt.min_score,
            opt.min_len))
        
        # Align each query along the subject an write result in a SAM file
        print ("Align queries against the subject sequence")
        i = 0
        for query in query_gen:
            
            # Find the alignment
            if opt.reverse:
                al, orient = find_best_align (ssw, query, float(opt.min_score), int(opt.min_len))
            else:
                al, orient = ssw.align(str(query.seq), float(opt.min_score), int(opt.min_len)), True
            
            # Try to extract the quality string if possible (fastq only)
            try:
                qual = SeqIO.QualityIO._get_sanger_quality_str(query)
            except ValueError:
                qual = "*"
            
            # Values if valid match found
            if al:
                flag = 0 if orient else 16
                rname = subject.id
                startpos = al.ref_begin
                cigar = al.cigar_string
                score = 0#al.score
            
            # Values if no valid match found
            else:
                flag = 4
                rname = '*'
                startpos = 0
                cigar = '*'
                score = 0

            # Line written in case of match
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                query.id,       # Read name
                flag,           # Flag 0 = mapped forward 16 = mapped forward reverse
                rname,          # Reference sequence
                startpos,       # Start position on reference
                score,          # MAPQ score = SSW score instead
                cigar,          # Cigar string
                '*',0,0,        # Pair end reads informations
                str(query.seq), # Sequence of read
                qual))          # Quality string of read (only for fastq query)
        
            # Progress bar 
            i+=1
            if i % 10 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
            if i % 1000 == 0:
                sys.stdout.write("\n")
                sys.stdout.flush()
                
        print ("\n{} Sequences processed in {}s".format(i, round(time()-start, 2)))

#~~~~~~~HELPER FUNCTIONS~~~~~~~#

def find_best_align (ssw, query, min_score, min_len):
    
    # Align reverse and forward query
    forward_al = ssw.align(str(query.seq), min_score, min_len)
    reverse_al = ssw.align(str(query.seq.reverse_complement()), min_score, min_len)
    
    # Decision tree to return the best aligned sequence taking into acount the absence of result
    # by ssw_wrap in case of score filtering
    
    if not forward_al:
        if not reverse_al:
            return (None, None)
        else:
            return (reverse_al, False)
    
    else:
        if not reverse_al:
            return (forward_al, True)
        else:
            if forward_al.score >= reverse_al.score:
                return (forward_al, True)
            else:
                return (reverse_al, False)
    
def optparser():
    
    print("Parse command line options")
    # Usage and version strings
    program_name = "pyssw"
    program_version = 0.1
    version_string = "{}\t{}".format(program_name, program_version)
    usage_string = "{}.py -s subject.fasta -q fastq (or fasta) [Facultative options]".format(program_name)
    optparser = optparse.OptionParser(usage = usage_string, version = version_string)

    # Define optparser options
    hstr = "Path of the fasta file containing the subject genome sequence [REQUIRED] "
    optparser.add_option( '-s', '--subject', dest="subject", help=hstr)
    hstr = "Path of the fastq or fasta file containing the short read to be aligned [REQUIRED]"
    optparser.add_option( '-q', '--query', dest="query", help=hstr)
    hstr = "positive integer for weight match in genome sequence alignment. [default: 2]"
    optparser.add_option( '-m', '--match', dest="match",default=2, help=hstr)
    hstr = "positive integer. The negative value will be used as weight mismatch in genome sequence alignment. [default: 2]"
    optparser.add_option( '-x', '--mismatch', dest="mismatch", default=2, help=hstr)
    hstr = "positive integer. The negative value will be used as weight for the gap opening. [default: 3]"
    optparser.add_option( '-o', '--gap_open', dest="gap_open", default=3, help=hstr)
    hstr = "positive integer. The negative value will be used as weight for the gap opening. [default: 1]"
    optparser.add_option( '-e', '--gap_extend', dest="gap_extend", default=1, help=hstr)
    hstr = "integer. Consider alignments having a score <= as not aligned. [default: 0]"
    optparser.add_option( '-f', '--min_score', dest="min_score", default=0, help=hstr)
    hstr = "integer. Consider alignments having a length <= as not aligned. [default: 0]"
    optparser.add_option( '-l', '--min_len', dest="min_len", default=0, help=hstr)
    hstr = "Flag. Align query in forward and reverse orientation and choose the best alignment. [default option]"
    optparser.add_option( '-r', '--reverse', dest="reverse", action="store_true", default=True, help=hstr)
    # Parse arg and return a dictionnary_like object of options
    opt, args = optparser.parse_args()
    
    if not opt.subject:
        print ("\nERROR: a subject fasta file has to be provided (-s option)\n")
        optparser.print_help()
        sys.exit()
    
    if not opt.query:
        print ("\nERROR: a query fasta or fastq file has to be provided (-q option)\n")
        optparser.print_help()
        sys.exit()
    
    return opt

#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#

if __name__ == '__main__':
    
    # try to import Third party and local packages   
    try:
        find_module('Bio')
        from Bio import SeqIO
    except ImportError:
        print ("ERROR: Please install Biopython package")
        sys.exit()
    
    try:
        find_module('ssw_wrap')
        from ssw_wrap import Aligner
    except ImportError:
        print ("ERROR: Please place ssw_wrap in the current directory or add its dir to python path")
        sys.exit()

    # Parse command line arguments
    opt = optparser()
    # Run the main function
    main(opt)

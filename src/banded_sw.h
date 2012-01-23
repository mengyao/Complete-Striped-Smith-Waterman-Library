/*
 *  banded_sw.h
 *  A banded Smith-Waterman algorithm to traceback the alignment path
 *  Created by Mengyao Zhao on 01/10/12.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 01/13/12.
 *
 */

/*! @function	A banded Smith-Waterman algorithm to traceback the alignment path 
 *  @parameter	read	sequence
 *  @return		alignment represented by a simple cigar string, regexp: ([1-9]+[MID])+
 */
char* banded_sw (const char* ref, 
				 	const char* read, 
				 	int32_t refLen, 
				 	int32_t readLen,
				 	uint32_t weight_match,    /* will be used as + */
				 	uint32_t weight_mismatch, /* will be used as - */
				 	uint32_t weight_insertB,  /* will be used as - */
				 	uint32_t weight_insertE,  /* will be used as - */
				 	uint32_t weight_deletB,   /* will be used as - */
				 	uint32_t weight_deletE,   /* will be used as - */
				 	int32_t band_width,
					int8_t* nt_table,
				 	int8_t* mat,			   /* pointer to the weight matrix */
				 	int32_t n);	

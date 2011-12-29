/*
 *  ssw.h
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 12/29/11.
 *	New features: Changed the data type of the score component of the alignment_end struct to unsigned char to avoid overflow.
 *
 */

#include <stdio.h>
#include <string.h>
#include <emmintrin.h>

/*! @typedef	struct of the alignment results
 *  @field score	the alignment score
 *	@field	ref	1-based position in the reference
 */
typedef struct {
	unsigned char score;
	int32_t ref;	/* 1-based position */
} alignment_end;

/*! @function	Generate query profile rearrange query sequence & calculate the weight of match/mismatch. 
 *  @parameter	read	sequence
 *  @parameter	weight_match	score for a pair of matched reference and read bases
 *  @parameter	weight_mismatch	score (absolute value) for a pair of mismached reference and read bases
 *	@parameter	bias	a number used to expend the max capacity of the values in the scoreing matrix; suggest to set to 4
 *  @return		pointer to the query profile 
 */
__m128i* queryProfile_constructor (const char* read,
								   unsigned char weight_match,    /* will be used as + */
								   unsigned char weight_mismatch, /* will be used as - */
								   unsigned char bias);

/*! @function	Transform the reference sequence to a number sequence. 
 *	@parameter	ref	reference sequence
 *	@parameter	refLen	reference length
 *	@return		reference sequence represented by numbers
 */
int32_t* ref_nt2num (const char* ref, int32_t refLen);


/*! @function	striped Smith-Waterman
 *  			Record the highest score of each reference position. 
 *  			Find the ending position of the optimal and sub-optimal alignment.
 *  @parameter	ref	reference sequence represented by numbers; can be generated using funceion ref_nt2num
 *  @parameter	refLen	reference length
 *  @parameter	readLen	read length
 *  @parameter	weight_insertB	score (absolute value) for opening a insertion 
 *  @parameter	weight_insertE	score (absolute value) for extending a insertion 
 *  @parameter	weight_deletB	score (absolute value) for opening a deletion 
 *  @parameter	weight_deletE	score (absolute value) for extending a deletion 
 *  @parameter	vProfile	pointer to the query profile
 *  @parameter	end_seg	a return value of 0-based segment number of alignment ending position in read; suggest to set to 0
 *	@parameter	bias	a number used to expend the max capacity of the values in the scoreing matrix; suggest to set to 4
 *  @return		a pointer to the array of structure alignment_end; the optimal (1st member of the array) and 
 *				sub-optimal (2nd member of the array) alignment score and ending positions
 */
alignment_end* smith_waterman_sse2 (const char* ref,
									int32_t refLen,
								    int32_t readLen, 
								    unsigned char weight_insertB, /* will be used as - */
								    unsigned char weight_insertE, /* will be used as - */
								    unsigned char weight_deletB,  /* will be used as - */
								    unsigned char weight_deletE,  /* will be used as - */
								    __m128i* vProfile,
								    int32_t* end_seg,        /* 0-based segment number of ending   
																	 alignment; The return value is  
																	 meaningful only when  
																     the return value != 0.
																   */	
	 							    unsigned char bias);	

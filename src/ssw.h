/*
 *  ssw.h
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 01/12/12.
 *	New features: Weight matrix is extracted.
 *
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <emmintrin.h>

/*! @typedef	struct of the alignment results
 *  @field score	the alignment score
 *	@field	ref	1-based position in the reference
 */
typedef struct {
	uint8_t score;
	int32_t ref;	/* 0-based position */
	int32_t read;	/* alignment ending position on read, 0-based */
} alignment_end;

/*! @function	Generate query profile rearrange query sequence & calculate the weight of match/mismatch. 
 *  @parameter	read	sequence
 *  @parameter	weight_match	score for a pair of matched reference and read bases
 *  @parameter	weight_mismatch	score (absolute value) for a pair of mismached reference and read bases
 *	@parameter	bias	a number used to expend the max capacity of the values in the scoreing matrix; suggest to set to 4
 *  @return		pointer to the query profile 
 */
__m128i* queryProfile_constructor (const char* read,
								   int8_t* mat,
								   int32_t n,	/* the edge length of the squre matrix mat */
								   uint8_t bias);

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
 *	@parameter	bias	a number used to expend the max capacity of the values in the scoreing matrix; suggest to set to 4
 *  @return		a pointer to the array of structure alignment_end; the optimal (1st member of the array) and 
 *				sub-optimal (2nd member of the array) alignment score and ending positions
 */
alignment_end* smith_waterman_sse2 (const char* ref,
									int32_t refLen,
								    int32_t readLen, 
								    uint8_t weight_insertB, /* will be used as - */
								    uint8_t weight_insertE, /* will be used as - */
								    uint8_t weight_deletB,  /* will be used as - */
								    uint8_t weight_deletE,  /* will be used as - */
								    __m128i* vProfile,
									uint8_t terminate,	/* the best alignment score: used to terminate 
														   the matrix calculation when locating the 
														   alignment beginning point. If this score 
														   is set to 0, it will not be used */
	 							    uint8_t bias);	

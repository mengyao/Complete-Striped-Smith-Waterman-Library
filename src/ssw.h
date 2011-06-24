/*
 *  ssw.cpp
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.3
 *	Revised by Mengyao Zhao on 4/5/11.
 *	New features: Record the highest score of each reference position. 
 *	Geve out the most possible 2nd distinguished best alignment score as well as the best alignment score and
 *	its ending position. 
 *
 */

#include <stdio.h>
#include <string.h>
#include <emmintrin.h>

/* struct of the alignment result */
typedef struct {
	char score;
	int32_t ref;	/* 1-based position */
} alignment_end;

/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch. */
__m128i* queryProfile_constructor (const char* read,
								   unsigned char weight_match,    /* will be used as + */
								   unsigned char weight_mismatch, /* will be used as - */
								   unsigned char bias);

/* Transform the reference sequence to a number sequence. */
int32_t* ref_nt2num (const char* ref, int32_t refLen);

/* striped Smith-Waterman
   Record the highest score of each reference position. 
   Find the ending position of the best alignment. 
   Beginning of gap and extention of gap are different. 
   wight_match > 0, all other weights < 0.
 */
alignment_end* smith_waterman_sse2 (const int32_t* ref,
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

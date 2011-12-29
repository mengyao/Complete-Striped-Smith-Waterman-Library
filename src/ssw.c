/*
 *  ssw.c
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 12/29/11.
 *	New features: Simplified the weighting matrix. 
 *
 */

#include <emmintrin.h>
#include "ssw.h"

/* This table is used to transform nucleotide letters into numbers. */
unsigned char nt_table[128] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 
};

/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch. */
__m128i* queryProfile_constructor (const char* read,
								 unsigned char weight_match,	/* will be used as + */
								 unsigned char weight_mismatch, /* will be used as - */
								 unsigned char bias) { 
					
	int32_t readLen = strlen(read);
	int32_t
	segLen = (readLen + 15) / 16; /* Split the 128 bit register into 16 pieces. 
								     Each piece is 8 bit. Split the read into 16 segments. 
								     Calculat 16 segments in parallel.
								   */
	__m128i* vProfile = (__m128i*)calloc(15 * segLen, sizeof(__m128i));
	int8_t* t = (int8_t*)vProfile;
	int32_t nt, i, j, k;
	int32_t segNum;
	int8_t mat[25];

	// initialize scoring matrix
	for (i = k = 0; i < 5; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? weight_match : -weight_mismatch;
		mat[k++] = 0; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = 0;
	
	/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch */
	for (nt = 0; nt < 5; nt ++) {
		for (i = 0; i < segLen; i ++) {
			j = i; 
			for (segNum = 0; segNum < 16 ; segNum ++) {
				*t++ = j>= readLen ? 0 : mat[nt * 5 + nt_table[(int)read[j]]] + bias;
				j += segLen;
			}
		}
	}

	return vProfile;
}

// Trace the max score from the max score vector.
unsigned char findMax (__m128i vMaxScore, 
					   int32_t* pos_vector) { // pos_vector should be originally 0.
	
	unsigned char max = 0;
	int32_t two[8];
	two[0] = _mm_extract_epi16(vMaxScore, 0);
	two[1] = _mm_extract_epi16(vMaxScore, 1);
	two[2] = _mm_extract_epi16(vMaxScore, 2);
	two[3] = _mm_extract_epi16(vMaxScore, 3);
	two[4] = _mm_extract_epi16(vMaxScore, 4);
	two[5] = _mm_extract_epi16(vMaxScore, 5);
	two[6] = _mm_extract_epi16(vMaxScore, 6);
	two[7] = _mm_extract_epi16(vMaxScore, 7);
	int32_t i;
	for (i = 0; i < 8; i ++) {
		int32_t low = two[i] & 0x00ff;
		if (low > max) {
			*pos_vector = i * 2;
			max = low;
		}
		int32_t high = two[i] & 0xff00;
		high = high >> 8;
		if (high > max) {
			*pos_vector = i * 2 + 1;
			max = high;
		}
	}
	return max;
}

/* Striped Smith-Waterman
   Record the highest score of each reference position. 
   Return the alignment score and ending position of the best alignment, 2nd best alignment, etc. 
   Gap begin and gap extention are different. 
   wight_match > 0, all other weights < 0.
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
	 							    unsigned char bias) {         /* Shift 0 point to a positive value. */
	
	*end_seg = 0; /* Initialize as aligned at num. 0 segment. */
	unsigned char max = 0;		                     /* the max alignment score */
	int32_t end_read = 0;
	int32_t end_ref = 0; /* 1_based best alignment ending point; Initialized as isn't aligned - 0. */
	int32_t segLen = (readLen + 15) / 16; /* number of segment */
	
	/* array to record the largest score of each reference position */
	char* maxColumn = (char*) calloc(refLen, 1); 
	
	/* array to record the alignment read ending position of the largest score of each reference position */
	int32_t* end_read_column = (int32_t*) calloc(refLen, sizeof(int32_t));
	
	/* Define 16 byte 0 vector. */
	__m128i vZero = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0); 
	vZero = _mm_xor_si128(vZero, vZero);

	__m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
	int32_t i;
	int32_t j;
	for (i = 0; i < segLen; i ++) {
		pvHStore[i] = vZero;
		pvE[i] = vZero;
	}
	
	/* 16 byte insertion begin vector */
	__m128i vInserB = _mm_set_epi8(weight_insertB, weight_insertB, weight_insertB, weight_insertB, 
								   weight_insertB, weight_insertB, weight_insertB, weight_insertB,
								   weight_insertB, weight_insertB, weight_insertB, weight_insertB, 
								   weight_insertB, weight_insertB, weight_insertB, weight_insertB);
	
	/* 16 byte insertion extention vector */
	__m128i vInserE = _mm_set_epi8(weight_insertE, weight_insertE, weight_insertE, weight_insertE, 
								   weight_insertE, weight_insertE, weight_insertE, weight_insertE, 
								   weight_insertE, weight_insertE, weight_insertE, weight_insertE, 
								   weight_insertE, weight_insertE, weight_insertE, weight_insertE);
	
	/* 16 byte deletion begin vector */
	__m128i vDeletB = _mm_set_epi8(weight_deletB, weight_deletB, weight_deletB, weight_deletB, 
								   weight_deletB, weight_deletB, weight_deletB, weight_deletB, 
								   weight_deletB, weight_deletB, weight_deletB, weight_deletB, 
								   weight_deletB, weight_deletB, weight_deletB, weight_deletB);

	/* 16 byte deletion extention vector */
	__m128i vDeletE = _mm_set_epi8(weight_deletE, weight_deletE, weight_deletE, weight_deletE, 
								   weight_deletE, weight_deletE, weight_deletE, weight_deletE, 
								   weight_deletE, weight_deletE, weight_deletE, weight_deletE, 
								   weight_deletE, weight_deletE, weight_deletE, weight_deletE);

	/* 16 byte bias vector */
	__m128i vBias = _mm_set_epi8(bias, bias, bias, bias, bias, bias, bias, bias,  
								 bias, bias, bias, bias, bias, bias, bias, bias);

	__m128i vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
	__m128i vMaxMark = vZero; /* Trace the highest score till the previous column. */	
	__m128i vTemp;
	
	int32_t edge;
	
	/* outer loop to process the reference sequence */
	for (i = 0; i < refLen; i ++) {
		
		int32_t cmp;
		__m128i vF = vZero; /* Initialize F value to 0. 
							   Any errors to vH values will be corrected in the Lazy_F loop. 
							 */
		__m128i vH = pvHStore[segLen - 1];
		vH = _mm_slli_si128 (vH, 1); /* Shift the 128-bit value in vH left by 1 byte. */
		
		/* Swap the 2 H buffers. */
		__m128i* pv = pvHLoad;
		__m128i vMarkColumn;
		
		__m128i vMaxColumn = vZero; /* vMaxColumn is used to record the max values of column i. */
		
		__m128i* vP = vProfile + nt_table[(int)ref[i]] * segLen; /* Right part of the vProfile */
		pvHLoad = pvHStore;
		pvHStore = pv;
		
		/* inner loop to process the query sequence */
		for (j = 0; j < segLen; j ++) {
			vH = _mm_adds_epu8(vH, vP[j]);
			vH = _mm_subs_epu8(vH, vBias); /* vH will be always > 0 */

			/* Get max from vH, vE and vF. */
			vH = _mm_max_epu8(vH, pvE[j]);
			vH = _mm_max_epu8(vH, vF);
			
			/* Update highest score encountered this far. */
			vMaxScore = _mm_max_epu8(vMaxScore, vH);
			vMaxColumn = _mm_max_epu8(vMaxColumn, vH);
			
			/* Save vH values. */
			pvHStore[j] = vH;

			/* Update vE value. */
			vH = _mm_subs_epu8(vH, vInserB); /* saturation arithmetic, result >= 0 */
			pvE[j] = _mm_subs_epu8(pvE[j], vInserE);
			pvE[j] = _mm_max_epu8(pvE[j], vH);
			
			/* Update vF value. */
			vF = _mm_subs_epu8(vF, vDeletE);
			vF = _mm_max_epu8(vF, vH);
			
			/* Load the next vH. */
			vH = pvHLoad[j];
		}
		
		/* Lazy_F loop
		   The computed vF value is for the given column. 
	       Since we are at the end, we need to shift the vF value over to the next column.
		 */
		vF = _mm_slli_si128 (vF, 1);
		
		/* Correct the vH values until there are no element in vF that could influence the vH values. */
		j = 0;
		vTemp = _mm_subs_epu8(pvHStore[j], vDeletB);
		vTemp = _mm_subs_epu8(vF, vTemp);
		vTemp = _mm_cmpeq_epi8(vTemp, vZero); /* result >= 0 */
		cmp = _mm_movemask_epi8(vTemp);
		
		while (cmp != 0xffff) {
			pvHStore[j] = _mm_max_epu8(pvHStore[j], vF);
			
			/* Update highest score incase the new vH value would change it. (New line I added!) */
			vMaxScore = _mm_max_epu8(vMaxScore, pvHStore[j]);
			vMaxColumn = _mm_max_epu8(vMaxColumn, pvHStore[j]);
			
			/* Update vE incase the new vH value would change it. */
			vH = _mm_subs_epu8(pvHStore[j], vInserB);
			pvE[j] = _mm_max_epu8(pvE[j], vH);
			
			/* Update vF value. */
			vF = _mm_subs_epu8(vF, vDeletE);
			
			j ++;
			if (j >= segLen) {
				j = 0;
				vF = _mm_slli_si128 (vF, 1);
			}
			vTemp = _mm_subs_epu8(pvHStore[j], vDeletB);
			vTemp = _mm_subs_epu8(vF, vTemp);
			vTemp = _mm_cmpeq_epi8(vTemp, vZero); /* result >= 0 */
			cmp = _mm_movemask_epi8(vTemp);
		}		
		/* end of Lazy-F loop */
		
		/* loop for tracing the ending point of the best alignment */
		vTemp = _mm_max_epu8(vMaxScore, vMaxMark);
		vTemp = _mm_cmpeq_epi8(vMaxMark, vTemp);
		cmp = _mm_movemask_epi8(vTemp);
		
		j = 0; 
		
		while (cmp != 0xffff) {
			vMaxMark = _mm_max_epu8(vMaxMark, pvHStore[j]);
			end_ref = i + 1; /* Adjust to 1-based position. */
			int32_t p = 0;
			unsigned char temp = findMax(vMaxMark, &p);
			
			if (temp > max) {
				*end_seg = j;
				end_read = p * segLen + j + 1;  
				max = temp;
			}
			
			vTemp = _mm_cmpeq_epi8(vMaxMark, vMaxScore);
			cmp = _mm_movemask_epi8(vTemp);
			j ++;
			if (j >= segLen) {
				j = 0;
			}
		}
		
		/* loop to record the max score of each column as well as the alignment ending position on the read */
		j = 0;
		maxColumn[i] = 0;
		vMarkColumn = pvHStore[j]; /* Record the cumulated max segment values. */
		cmp = _mm_movemask_epi8(vMarkColumn);
		while (cmp != 0xffff) {
			vMarkColumn = _mm_max_epu8(vMarkColumn, pvHStore[j]);
			
			int32_t p = 0;
			unsigned char temp = findMax(vMarkColumn, &p);
			if (temp > maxColumn[i]) {
				end_read_column[i] = p * segLen + j + 1;  
				maxColumn[i] = temp;
			}
			vTemp = _mm_cmpeq_epi8(vMarkColumn, vMaxColumn);
			cmp = _mm_movemask_epi8(vTemp);
			j ++;
			if (j >= segLen) {
				j = 0;
			}
		}	
		
	} 	
	
	/* Find the most possible 2nd best alignment. */
	alignment_end* bests = (alignment_end*) calloc(2, sizeof(alignment_end));
	bests[0].score = max;
	bests[0].ref = end_ref;
	
	bests[1].score = 0;
	bests[1].ref = 0;

	edge = (end_ref - readLen / 2 - 1) > 0 ? (end_ref - readLen / 2 - 1) : 0;
	for (i = 0; i < edge; i ++) {
		if (maxColumn[i] > bests[1].score) {
			bests[1].score = maxColumn[i];
			bests[1].ref = i + 1;
		}
	}
	edge = (end_ref + readLen / 2 + 1) > refLen ? refLen : (end_ref + readLen / 2 + 1);
	for (i = edge; i < refLen; i ++) {
		if (maxColumn[i] > bests[1].score) {
			bests[1].score = maxColumn[i];
			bests[1].ref = i + 1;
		}		
	}
	
	free(pvHStore);
	free(maxColumn);
	free(end_read_column);
	return bests;
}


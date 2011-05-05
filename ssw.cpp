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

using namespace std;
#include <emmintrin.h>
#include "ssw.h"

// for Weight matrix
#define A 0
#define C 1
#define G 2
#define T 3
#define K 4  // G or T
#define M 5  // A or C
#define R 6  // A or G
#define Y 7  // C or T
#define S 8  // A or T
#define B 9  // C or G or T
#define V 10 // A or C or G
#define H 11 // A or C or T
#define D 12 // A or G or T
#define N 13 // any
#define X 14 // any, X mask on Y chromosome

// Transform amino acid to numbers for Weight Matrix.
unsigned int amino2num (char amino) {
	unsigned int num;
	switch (amino) {
		case 'A':
		case 'a':
			num = 0;
			break;
		case 'C':
		case 'c':
			num = 1;
			break;
		case 'G':
		case 'g':
			num = 2;
			break;
		case 'T':
		case 't':
			num = 3;
			break;
		case 'K': // G or T
		case 'k':
			num = 4;
			break;
		case 'M': // A or C
		case 'm':
			num = 5;
			break;
		case 'R': // A or G
		case 'r':
			num = 6;
			break;
		case 'Y': // C or T
		case 'y':
			num = 7;
			break;
		case 'S': // A or T
		case 's':
			num = 8;
			break;
		case 'B': // C or G or T
		case 'b':
			num = 9;
			break;
		case 'V': // A or C or G
		case 'v':
			num = 10;
			break;
		case 'H': // A or C or T
		case 'h':
			num = 11;
			break;
		case 'D': // A or G or T
		case 'd':
			num = 12;
			break;
		case 'N': // any
		case 'n':
			num = 13;
			break;
		case 'X': // any, x mask on Y chromosome
		case 'x':
			num = 14;
			break;
		default:
			fprintf(stderr, "Wrong sequence. \n");
			break;
	}
	return num;
}

//Transform the reference sequence to a number sequence.
unsigned int* ref_amino2num (const char* ref, unsigned int refLen) {
	unsigned int* ref_num = (unsigned int *)calloc(refLen, sizeof(unsigned int));
	for (unsigned int i = 0; i < refLen; i ++) {
		ref_num[i] = amino2num(ref[i]);
	}
	return ref_num;
}

// Create scoring matrix W.
char** matrixScore_constructor (unsigned char weight_match, 
								unsigned char weight_mismatch) {
	
	char** W = (char**)calloc(15, sizeof(char*));
	for (unsigned int i = 0; i < 15; i ++) {
		W[i] = (char*) calloc(15, sizeof(char));
	}
	for (unsigned int i = 0; i < 15; i ++) {
		for (unsigned int j = 0; j < 15; j ++) {
			if (i == j) {
				W[i][j] = weight_match;
			} else {
				W[i][j] = - weight_mismatch;
			}
		}
	}
	
	W[G][K] = W[K][G] = weight_match; // K
	W[T][K] = W[K][T] = weight_match;
	
	W[A][M] = W[M][A] = weight_match; // M
	W[C][M] = W[M][C] = weight_match;
	
	W[A][R] = W[R][A] = weight_match; // R
	W[G][R] = W[R][G] = weight_match;
	
	W[C][Y] = W[Y][C] = weight_match; // Y
	W[T][Y] = W[Y][T] = weight_match;
	
	W[A][S] = W[S][A] = weight_match; // S
	W[T][S] = W[S][T] = weight_match;
	
	W[C][B] = W[B][C] = weight_match; // B
	W[G][B] = W[B][G] = weight_match;
	W[T][B] = W[B][T] = weight_match;
	
	W[A][V] = W[V][A] = weight_match; // V
	W[C][V] = W[V][C] = weight_match; 
	W[G][V] = W[V][G] = weight_match;
	
	W[A][H] = W[H][A] = weight_match; // H
	W[C][H] = W[H][C] = weight_match; 
	W[T][H] = W[H][T] = weight_match;
	
	W[A][D] = W[D][A] = weight_match; // D
	W[G][D] = W[D][G] = weight_match;
	W[T][D] = W[D][T] = weight_match;
	
	return W;
}

// Generate query profile：rearrange query sequence & calculate the weight of match/mismatch.
__m128i* queryProfile_constructor (const char* read,
								 unsigned char weight_match,	// will be used as +
								 unsigned char weight_mismatch, // will be used as -
								 unsigned char bias) { 
					
	char** W = matrixScore_constructor (weight_match, weight_mismatch); // Create scoring matrix.
	
	// Generate query profile：rearrange query sequence & calculate the weight of match/mismatch.
	unsigned int readLen = strlen(read);
	unsigned int
	segLen = (readLen + 15) / 16; // Split the 128 bit register into 16 pieces. 
								  // Each piece is 8 bit. Split the read into 16 segments. 
								  // Calculat 16 segments in parallel.
	__m128i* vProfile = (__m128i*)calloc(15 * segLen, sizeof(__m128i));
	
	int8_t* t = (int8_t*)vProfile;
	for (unsigned int amino = 0; amino < 15; amino ++) {
			for (unsigned int i = 0; i < segLen; i ++) {
			for (unsigned int j = i; j < 16; j += segLen) {
				*t++ = j>= readLen ? 0 : W[amino][amino2num(read[j])] + bias;
			}
		}
	}
	
	for (unsigned int i = 0; i < 15; i ++) {
		free(W[i]);
	}
	free(W);
	return vProfile;
}

// Striped Smith-Waterman
// Record the highest score of each reference position. 
// Return the alignment score and ending position of the best alignment, 2nd best alignment, etc. 
// Gap begin and gap extention are different. 
// wight_match > 0, all other weights < 0. 
alignment_end* smith_waterman_sse2 (const unsigned int* ref_num,
									unsigned int refLen,
								    unsigned int readLen, 
								    unsigned char weight_insertB, // will be used as -
								    unsigned char weight_insertE, // will be used as -
								    unsigned char weight_deletB,  // will be used as -
								    unsigned char weight_deletE,  // will be used as -
								    __m128i* vProfile,
								    unsigned int* end_seg,        // 0-based segment number of ending  
								                                  // alignment; The return value is  
																  // meaningful only when  
									  							  // the return value != 0.
	 							    unsigned char bias) {         // Shift 0 point to a positive value.
	
	*end_seg = 0; // Initialize as aligned at num. 0 segment.
	unsigned char max = 0;		                     // the max alignment score
	unsigned int end_ref = 0; // 1_based best alignment ending point; Initialized as isn't aligned - 0.
	unsigned int segLen = (readLen + 15) / 16; // number of segment
	
	// array to record the largest score of each reference position
	char* maxColumn = (char*) calloc(refLen, 1); 
	
	// array to record the alignment read ending position of the largest score of each reference position 
	unsigned int* end_read_column = (unsigned int*) calloc(refLen, sizeof(unsigned int));
	
	// Define 16 byte 0 vector.
	__m128i vZero; 
	vZero = _mm_xor_si128(vZero, vZero);

	__m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
	for (unsigned int i = 0; i < segLen; i ++) {
		pvHStore[i] = vZero;
		pvE[i] = vZero;
	}
	
	// 16 byte insertion begin vector
	__m128i vInserB = _mm_set_epi8(weight_insertB, weight_insertB, weight_insertB, weight_insertB, 
								   weight_insertB, weight_insertB, weight_insertB, weight_insertB,
								   weight_insertB, weight_insertB, weight_insertB, weight_insertB, 
								   weight_insertB, weight_insertB, weight_insertB, weight_insertB);
	
	// 16 byte insertion extention vector
	__m128i vInserE = _mm_set_epi8(weight_insertE, weight_insertE, weight_insertE, weight_insertE, 
								   weight_insertE, weight_insertE, weight_insertE, weight_insertE, 
								   weight_insertE, weight_insertE, weight_insertE, weight_insertE, 
								   weight_insertE, weight_insertE, weight_insertE, weight_insertE);
	
	// 16 byte deletion begin vector
	__m128i vDeletB = _mm_set_epi8(weight_deletB, weight_deletB, weight_deletB, weight_deletB, 
								   weight_deletB, weight_deletB, weight_deletB, weight_deletB, 
								   weight_deletB, weight_deletB, weight_deletB, weight_deletB, 
								   weight_deletB, weight_deletB, weight_deletB, weight_deletB);

	// 16 byte deletion extention vector
	__m128i vDeletE = _mm_set_epi8(weight_deletE, weight_deletE, weight_deletE, weight_deletE, 
								   weight_deletE, weight_deletE, weight_deletE, weight_deletE, 
								   weight_deletE, weight_deletE, weight_deletE, weight_deletE, 
								   weight_deletE, weight_deletE, weight_deletE, weight_deletE);

	// 16 byte bias vector
	__m128i vBias = _mm_set_epi8(bias, bias, bias, bias, bias, bias, bias, bias,  
								 bias, bias, bias, bias, bias, bias, bias, bias);

	__m128i vMaxScore = vZero; // Trace the highest score of the whole SW matrix.
	__m128i vMaxMark = vZero; // Trace the highest score till the previous column.	
	__m128i vTemp;
	
	// outer loop to process the reference sequence
	for (unsigned int i = 0; i < refLen; i ++) {
		
		__m128i vF = vZero; // Initialize F value to 0. 
							// Any errors to vH values will be corrected in the Lazy_F loop. 
		__m128i vH = pvHStore[segLen - 1];
		vH = _mm_slli_si128 (vH, 1); // Shift the 128-bit value in vH left by 1 byte.
		
		// Swap the 2 H buffers.
		__m128i* pv = pvHLoad;
		pvHLoad = pvHStore;
		pvHStore = pv;
		
		__m128i vMaxColumn = vZero; // vMaxColumn is used to record the max values of column i.
		
		__m128i* vP = vProfile + ref_num[i] * segLen; //Right part of the vProfile
		
		// inner loop to process the query sequence
		for (unsigned int j = 0; j < segLen; j ++) {
			vH = _mm_adds_epu8(vH, vP[j]);
			vH = _mm_subs_epu8(vH, vBias); // vH will be always > 0

			// Get max from vH, vE and vF.
			vH = _mm_max_epu8(vH, pvE[j]);
			vH = _mm_max_epu8(vH, vF);
			
			// Update highest score encountered this far.
			vMaxScore = _mm_max_epu8(vMaxScore, vH);
			vMaxColumn = _mm_max_epu8(vMaxColumn, vH);
			
			// Save vH values.
			pvHStore[j] = vH;

			// Update vE value.
			vH = _mm_subs_epu8(vH, vInserB); // saturation arithmetic, result >= 0
			pvE[j] = _mm_subs_epu8(pvE[j], vInserE);
			pvE[j] = _mm_max_epu8(pvE[j], vH);
			
			// Update vF value.
			vF = _mm_subs_epu8(vF, vDeletE);
			vF = _mm_max_epu8(vF, vH);
			
			// Load the next vH.
			vH = pvHLoad[j];
		}
		
		// Lazy_F loop
		// The computed vF value is for the given column. 
		// Since we are at the end, we need to shift the vF value over to the next column.
		vF = _mm_slli_si128 (vF, 1);
		
		// Correct the vH values until there are no element in vF that could influence the vH values.
		unsigned int j = 0;
		vTemp = _mm_subs_epu8(pvHStore[j], vDeletB);
		vTemp = _mm_subs_epu8(vF, vTemp);
		vTemp = _mm_cmpeq_epi8(vTemp, vZero); // result >= 0
		unsigned int cmp = _mm_movemask_epi8(vTemp);
		
		while (cmp != 0xffff) {
			pvHStore[j] = _mm_max_epu8(pvHStore[j], vF);
			
			// Update highest score incase the new vH value would change it. (New line I added!)
			vMaxScore = _mm_max_epu8(vMaxScore, pvHStore[j]);
			vMaxColumn = _mm_max_epu8(vMaxColumn, pvHStore[j]);
			
			// Update vE incase the new vH value would change it.
			vH = _mm_subs_epu8(pvHStore[j], vInserB);
			pvE[j] = _mm_max_epu8(pvE[j], vH);
			
			// Update vF value.
			vF = _mm_subs_epu8(vF, vDeletE);
			
			j ++;
			if (j >= segLen) {
				j = 0;
				vF = _mm_slli_si128 (vF, 1);
			}
			vTemp = _mm_subs_epu8(pvHStore[j], vDeletB);
			vTemp = _mm_subs_epu8(vF, vTemp);
			vTemp = _mm_cmpeq_epi8(vTemp, vZero); // result >= 0
			cmp = _mm_movemask_epi8(vTemp);
		}		
		// end of Lazy-F loop
		
		// loop for tracing the ending point of the best alignment
		//vTemp = vMaxMark;
		vTemp = _mm_max_epu8(vMaxScore, vMaxMark);
		vTemp = _mm_cmpeq_epi8(vMaxMark, vTemp);
		cmp = _mm_movemask_epi8(vTemp);
		
		j = 0; 
		
		while (cmp != 0xffff) {
			vMaxMark = _mm_max_epu8(vMaxMark, pvHStore[j]);
			end_ref = i + 1; // Adjust to 1-based position.
			
			/* find largest score in the vMaxMark vector */
			__m128i vT = _mm_srli_si128 (vMaxMark, 8);
			__m128i vM = _mm_max_epi16 (vMaxMark, vT);
			vT = _mm_srli_si128 (vM, 4);
			vM = _mm_max_epi16 (vM, vT);
			vT = _mm_srli_si128 (vM, 2);
			vM = _mm_max_epi16 (vM, vT);
			vT = _mm_srli_si128 (vM, 1);
			vM = _mm_max_epi16 (vM, vT);
			unsigned int temp = _mm_extract_epi16 (vM, 0);
			
			if (temp > max) {
				*end_seg = j;
				max = temp;
			}
			
			vTemp = _mm_cmpeq_epi8(vMaxMark, vMaxScore);
			cmp = _mm_movemask_epi8(vTemp);
			j ++;
			if (j >= segLen) {
				j = 0;
			}
		}
		
		// loop to record the max score of each column as well as the alignment ending position on the read
		j = 0;
		maxColumn[i] = 0;
		__m128i vMarkColumn = pvHStore[j]; // Record the cumulated max segment values.
		cmp = _mm_movemask_epi8(vMarkColumn);
		while (cmp != 0xffff) {
			
			vMarkColumn = _mm_max_epu8(vMarkColumn, pvHStore[j]);
			
			/* find largest score in the vMarkColumn vector */
			__m128i vT = _mm_srli_si128 (vMarkColumn, 8);
			__m128i vM = _mm_max_epi16 (vMarkColumn, vT);
			vT = _mm_srli_si128 (vM, 4);
			vM = _mm_max_epi16 (vM, vT);
			vT = _mm_srli_si128 (vM, 2);
			vM = _mm_max_epi16 (vM, vT);
			vT = _mm_srli_si128 (vM, 1);
			vM = _mm_max_epi16 (vM, vT);
			unsigned int temp = _mm_extract_epi16 (vM, 0);
			
			if (temp > maxColumn[i]) {
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
	
	// Find the most possible 2nd best alignment.
	alignment_end* bests = (alignment_end*) calloc(2, sizeof(alignment_end));
	bests[0].score = max;
	bests[0].ref = end_ref;
	
	bests[1].score = 0;
	bests[1].ref = 0;

	unsigned int edge = (end_ref - readLen / 2 - 1) > 0 ? (end_ref - readLen / 2 - 1) : 0;
	for (unsigned int i = 0; i < edge; i ++) {
		if (maxColumn[i] > bests[1].score) {
			bests[1].score = maxColumn[i];
			bests[1].ref = i + 1;
		}
	}
	edge = (end_ref + readLen / 2 + 1) > refLen ? refLen : (end_ref + readLen / 2 + 1);
	for (unsigned int i = edge; i < refLen; i ++) {
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


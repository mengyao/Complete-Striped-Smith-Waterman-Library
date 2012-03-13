/*
 *  ssw.h
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 03/12/12.
 *	New features: This is the api file.
 *
 */

#ifndef SSW_H
#define SSW_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <emmintrin.h>

struct _profile;
typedef struct _profile s_profile;

// Positions are all 1-based.
typedef struct {
	uint16_t score1;	// best alignment score, 225: best alignment score is > 225
	uint16_t score2;	// sub-optimal alignment score
	int32_t ref_begin1;	// best alignment beginning position on reference, 0: none
	int32_t ref_end1;	// best alignment ending position on reference
	int32_t	read_begin1;	// best alignment beginning position on read, 0: none
	int32_t read_end1;	// best alignment ending position on read
	int32_t ref_end2;	// sub-optimal alignment ending position on reference
	uint32_t* cigar;	// best alignment cigar, the same as that in Bam format, high 28 bits: length, low 4 bits: M/I/D (0/1/2), 0: none
	int32_t cigarLen;	// length of the cigar string
} s_align;

#ifdef __cplusplus
extern "C" {
#endif	// __cplusplus

// @function	Create the ssw profile using the read sequence.
s_profile* ssw_init (const int8_t* read, const int32_t readLen, const int8_t* mat, const int32_t n, const int8_t score_size);

// @function	Release the memory alloced by function ssw_init.
void init_destroy (s_profile* p);

// @function	ssw alignment.
s_align* ssw_align (const s_profile* prof, 
					const int8_t* ref, 
					const int32_t refLen, 
					const uint8_t weight_gapO, 
					const uint8_t weight_gapE, 
					const uint8_t flag,	//  (from high to low) bit 6: return the best alignment beginning position; 7: if max score >= filter, return cigar; 8: always return cigar
					const uint16_t filter);

// @function	Release the memory alloced by function ssw_align.
void align_destroy (s_align* a);

#ifdef __cplusplus
}
#endif	// __cplusplus

#endif	// SSW_H

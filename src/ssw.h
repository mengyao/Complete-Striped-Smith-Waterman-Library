/*
 *  ssw.h
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 03/26/12.
 *	New features: This is the api file.
 *
 */

#ifndef SSW_H
#define SSW_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <emmintrin.h>

/*!	@typedef	structure of the return values of the function ssw_init; contains the quiry profile	*/
struct _profile;
typedef struct _profile s_profile;

/*!	@typedef	structure of the return values of the function ssw_align
	@field	score1	the best alignment score; score1 = 255 when the best alignment score is >= 255
	@field	score2	sub-optimal alignment score
	@field	ref_begin1	best alignment beginning position on reference;	ref_begin1 = 0 when the best alignment beginning position 
						is not available
	@field	ref_end1	best alignment ending position on reference
	@field	read_begin1	best alignment beginning position on read; read_begin1 = 0 when the best alignment beginning position is 
						not available
	@field	read_end1	best alignment ending position on read
	@field	read_end2	sub-optimal alignment ending position on read
	@field	cigar	best alignment cigar; stored the same as that in BAM format, high 28 bits: length, low 4 bits: M/I/D (0/1/2); 
					cigar = 0 when the best alignment path is not available
	@field	cigarLen	length of the cigar string; cigarLen = 0 when the best alignment path is not available
	@note	The fields ref_begin1, ref_end1, read_begin1 read_end1 and read_end2 all have 1-based coordinate.
*/
typedef struct {
	uint16_t score1;	
	uint16_t score2;	
	int32_t ref_begin1;	
	int32_t ref_end1;	
	int32_t	read_begin1;	
	int32_t read_end1;	
	int32_t ref_end2;
	uint32_t* cigar;	
	int32_t cigarLen;	
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
					int32_t refLen, 
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

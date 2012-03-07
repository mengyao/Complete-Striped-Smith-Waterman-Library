/*
 *  ssw.h
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 03/07/12.
 *	New features: This is the api file.
 *
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <emmintrin.h>

struct _profile;
typedef struct _profile s_profile;

// Positions are all 1-based.
typedef struct {
	int16_t score1;	// best alignment score, 225: best alignment score is > 225
	int16_t score2;	// sub-optimal alignment score
	int32_t ref_begin1;	// best alignment beginning position on reference, 0: none
	int32_t ref_end1;	// best alignment ending position on reference
	int32_t	read_begin1;	// best alignment beginning position on read, 0: none
	int32_t read_end1;	// best alignment ending position on read
	int32_t ref_end2;	// sub-optimal alignment ending position on reference
	char* cigar;	// best alignment cigar, 0: none
} s_align;
/*
typedef struct {
	align* al;
	char* ref_name;
	char* ref_seq;
	char* read_name;
	char* read_seq;	// strand == 0: original read; strand == 1: reverse complement read
	int8_t strand;	// 0: forward aligned ; 1: reverse complement aligned
	int8_t sam;	// 0: Blast like output; 1: Sam format output
} write_param;
*/
// @function	Create the ssw profile using the read sequence.
s_profile* ssw_init (int8_t* read, int32_t readLen, int8_t* mat, int32_t n, int8_t score_size);

// @function	Release the memory alloced by function ssw_init.
void init_destroy (s_profile* p);

// @function	ssw alignment.
s_align* ssw_align (s_profile* prof, int8_t* ref, int32_t refLen, uint8_t weight_gapO, uint8_t weight_gapE, int8_t begin, int8_t align);

// @function	Release the memory alloced by function ssw_align.
void align_destroy (s_align* c);

//void ssw_write (write_param* w);

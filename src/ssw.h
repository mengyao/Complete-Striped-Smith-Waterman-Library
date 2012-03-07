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

typedef struct {
	int8_t* read;	// read sequence represented in numbers
	int8_t* mat;	// substitution matrix correspounds to the read sequence
	int8_t score_size;	// 0: best alignment score will be < 225; 1: > 225; 2: can be either
	int32_t readLen;	// read length
	int32_t n;	// column number of the subsitution matrix mat
} init_param;

struct _profile;
typedef struct _profile profile;

typedef struct {
	profile* prof;
	int8_t* ref;	// reference sequence represented in numbers
	int32_t refLen;	// reference length
	uint8_t weight_gapO; /* will be used as - */
	uint8_t weight_gapE; /* will be used as - */
//	uint8_t weight_deletB;  /* will be used as - */
//	uint8_t weight_deletE;  /* will be used as - */
	int8_t begin;	// 1: the best alignment beginning position is needed; 0: otherwise
	int8_t align;	// 1: the best alignment path (cigar) is needed; 0: otherwise
} align_param;

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
} align;

typedef struct {
	align* al;
	char* ref_name;
	char* ref_seq;
	char* read_name;
	char* read_seq;	// strand == 0: original read; strand == 1: reverse complement read
	int8_t strand;	// 0: forward aligned ; 1: reverse complement aligned
	int8_t sam;	// 0: Blast like output; 1: Sam format output
} write_param;

// @function	Create the ssw profile using the read sequence.
profile* ssw_init (init_param* init);

// @function	Release the memory alloced by function ssw_init.
void init_destroy (profile* p);

// @function	ssw alignment.
align* ssw_align (align_param* a);

// @function	Release the memory alloced by function ssw_align.
void align_destroy (align* c);

//void ssw_write (write_param* w);

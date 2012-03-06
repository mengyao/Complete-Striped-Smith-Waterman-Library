/*
 *  ssw.h
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 03/05/12.
 *	New features: This is the api file.
 *
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <emmintrin.h>

typedef struct {
	int8_t* read;
	int8_t* rc_read;	// 0: reverse complement alignment will not be done; vise versa
	int8_t* mat;
	int8_t score_size;	// 0: best alignment score will be < 225; 1: > 225; 2: can be either
	int32_t readLen;
	int32_t n;
} init_param;

struct _profile;
typedef struct _profile profile;

typedef struct {
	profile* prof;
	int8_t* ref;
	int32_t refLen;
	uint8_t weight_insertB; /* will be used as - */
	uint8_t weight_insertE; /* will be used as - */
	uint8_t weight_deletB;  /* will be used as - */
	uint8_t weight_deletE;  /* will be used as - */
	int8_t begin;	// 1: the best alignment beginning position is needed; 0: otherwise
	int8_t align;	// 1: the best alignment path (cigar) is needed; 0: otherwise
} align_param;

// Positions are all 1-based.
typedef struct {
	int8_t strand;	// 0: forward aligned; 1: reverse complement aligned 
	int16_t score1;	// best alignment score, 225: 
	int16_t score2;	// sub-optimal alignment score
	int32_t ref_begin1;	// 0: none
	int32_t ref_end1;
	int32_t	read_begin1;	// 0: none
	int32_t read_end1;
	int32_t ref_end2;
	char* cigar;	// best alignment cigar, 0: none
} align;

profile* ssw_init (init_param* init);

void init_destroy (profile* p);

align* ssw_align (align_param* a);

void align_destroy (align* c);

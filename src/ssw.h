/*
 *  ssw.h
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 03/03/12.
 *	New features: This is the api file.
 *
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <emmintrin.h>

typedef struct {
	const char* read;
	int8_t* mat;
	int8_t score_size;	// 0: best alignment score will be < 225; 1: > 225; 2: can be either
	int8_t type;	// 0: genome sequence; 1: protein sequence
	int8_t reverse;	// 1: reverse complement alignment will also be done; 0: otherwise
} init_param;

struct _profile;
typedef struct _profile profile;

typedef struct {
	profile* prof;
	const char* ref;
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
	const char* read;	// if strand == 0: original read seq; else reverse complement read seq
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

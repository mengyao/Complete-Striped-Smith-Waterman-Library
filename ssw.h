/*
 *  ssw.cpp
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.3
 *	Revised by Mengyao Zhao on 8/27/10.
 *	New features: Record the highest score of each reference position. 
 *	Geve out the most possible 2nd distinguished best alignment score as well as the best alignment score and
 *	its ending position. 
 *
 */

using namespace std;
#include <vector>
#include <stdio.h>
#include <string.h>



// struct of the alignment result
typedef struct {
	char score;
	unsigned int ref;	// 1-based position
//	unsigned int read;	// 1-based position
} alignment_end;

// Generate query profile：rearrange query sequence & calculate the weight of match/mismatch.
char** queryProfile_constructor (const char* read,
								 unsigned char weight_match,    // will be used as +
								 unsigned char weight_mismatch, // will be used as -
								 unsigned char bias);

// Free the memory of queryProfile.
void queryProfile_destructor (char** queryProfile);

//Transform the reference sequence to a number sequence.
unsigned int* ref_amino2num (const char* ref, unsigned int refLen);

// striped Smith-Waterman
// Record the highest score of each reference position. 
// Find the ending position of the best alignment. 
// Beginning of gap and extention of gap are different. 
// wight_match > 0, all other weights < 0.
alignment_end* smith_waterman_sse2 (const unsigned int* ref,
									unsigned int refLen,
								    unsigned int readLen, 
								    unsigned char weight_insertB, // will be used as -
								    unsigned char weight_insertE, // will be used as -
								    unsigned char weight_deletB,  // will be used as -
								    unsigned char weight_deletE,  // will be used as -
								    char** queryProfile,
								    unsigned int* end_seg,        // 0-based segment number of ending  
									// alignment; The return value is  
									// meaningful only when  
									// the return value != 0.
	 							    unsigned char bias);	

// Reverse Smith-Waterman beginning from the end of the best alignment. 
// Alignment beginning position will be signed to end_ref & end_read when finish running.
void reverse_smith_waterman_sse2 (const char* ref, 
								  unsigned int segLen, 
								  unsigned char weight_insertB, // will be used as -
								  unsigned char weight_insertE, // will be used as -
								  unsigned char weight_deletB, // will be used as -
								  unsigned char weight_deletE, // will be used as -
								  char** queryProfile,
								  unsigned int end_ref, // 1-based position of the alignment ending. 0 is not aligned.
								  unsigned int* begin_ref, // 1-based position of the alignment begin
								  unsigned int* begin_read, // 1-based position of the alignment begin
								  unsigned int end_seg, // 0-based segment number of ending alignment
								  unsigned char max, // the best alignment score 
								  unsigned char bias); // Shift 0 point to a positive value.
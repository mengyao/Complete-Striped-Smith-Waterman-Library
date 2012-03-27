/*
 *  ssw.h
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 03/27/12.
 *	New features: This is the api file.
 *
 */

#ifndef SSW_H
#define SSW_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <emmintrin.h>

/*!	@typedef	structure of the return values of the function ssw_init; contains the query profile	*/
struct _profile;
typedef struct _profile s_profile;

/*!	@typedef	structure of the alignment result
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

/*!	@function	Create the query profile using the query sequence.
	@param	read	pointer to the query sequence; the query sequence is represented in numbers
	@param	readLen	length of the query sequence
	@param	mat	pointer to the substitution matrix
	@param	n	the number of elements in mat is n*n
	@param	score_size	estimated Smith-Waterman score; if your estimated best alignment score is surely < 255 or you would like t 
						stop the best alignment searching when its score reaches 255, please set 0; if your estimated best 
						alignment score >= 255, please set 1; if you don't know, please set 2 
	@return	pointer to the query profile structure
*/
s_profile* ssw_init (const int8_t* read, const int32_t readLen, const int8_t* mat, const int32_t n, const int8_t score_size);

/*!	@function	Release the memory allocated by function ssw_init.
	@param	p	pointer to the query profile structure	
*/
void init_destroy (s_profile* p);

// @function	ssw alignment.
/*!	@function	Do Striped Smith-Waterman alignment.
	@param	prof	pointer to the query profile structure
	@param	ref	pointer to the target sequence; the target sequence is represented in numbers
	@param	refLen	length of the target sequence
	@param	weight_gapO	the absolute value of gap open penalty  
	@param	weight_gapE	the absolute value of gap extension penalty
	@param	flag	bitwise FLAG; (from high to low) bit 6: when setted as 1, function ssw_align will return the best alignment 
					beginning position; bit 7: when setted as 1, if the best alignment score >= filter, (whatever bit 6 is setted)
					the function will return the best alignment beginning position and cigar; bit 8: when setted as 1, (whatever 
					bit 6 or 7 is setted) the function will always return the best alignment beginning position and cigar
	@param	filter	when bit 7 of flag is setted as 1 and bit 8 is setted as 0, filter will be used
	@return	pointer to the alignment result structure  	
*/
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

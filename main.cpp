/*
 *  ssw.cpp
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.3
 *	Revised by Mengyao Zhao on 8/27/10.
 *	New features: Record the highest 2*read_length scores among the highest score of each reference position. 
 *	Geve out the most possible 2nd distinguished best alignment score as well as the best alignment score and
 *	its ending position. 
 *
 */

using namespace std;

#include <stdlib.h>
#include <vector>
#include "ssw.h"

// For BandedSmithWaterman:
#include <iostream>
#include <string.h>
#include <string>

//Read reference and read sequences from the test file.
void readFasta (char* ref_read, char* ref, char* read) {
	FILE* Handle = fopen (ref_read, "r");
	if(Handle == NULL) {
		fprintf(stderr, "Error: Can't access the test fasta file with reference and read sequences: %s.\n", 
				ref_read);
	}
	unsigned int s = SIZE, i = 0, flag = 0; //flag == 1: reference sequence; flag == 2: read
	char* line = (char*) calloc(SIZE, 1);
	while (! feof (Handle)) {
		char letter = fgetc (Handle);
		if (letter == '\n') {
			line[i] = '\0';
			if (strstr(line, ">ref")) {
				flag = 1;
			} else if (strstr(line, ">read")){
				flag = 2;
			} else if (i != 0){
				if (flag == 1) {
					if (i == SIZE - 1) {
						ref = (char*) realloc (ref, i + 1);
					}
					strcpy (ref, line);
				} else if (flag == 2) {
					if (i == SIZE - 1) {
						read = (char*) realloc (read, i + 1);
					}
					strcpy (read, line);
				}
			}
			i = 0;
		} else {
			if (i == s - 1) {
				s *= 2;
				line = (char*) realloc (line, s);
			}
			line[i] = letter;
			i ++;
		}
	}
	fclose(Handle);
}

int main (int argc, char * const argv[]) {
	// validate argument count
	if (argc < 2) {
		fprintf (stderr, "USAGE: %s <ref_read.fa>\n", argv[0]);
		exit (1);
	}
	char* ref_read = argv[1];
	char* ref = (char*) calloc(SIZE, 1); 
	char* read = (char*) calloc(SIZE, 1);
	readFasta(ref_read, ref, read);
	
	char** queryProfile = queryProfile_constructor(read, 2, 1, 4);
	unsigned int end_seg = 0;
	unsigned int readLen = strlen(read);
	alignment_end* bests = smith_waterman_sse2(ref, readLen, 2, 1, 2, 1, queryProfile, &end_seg, 4);
/*	reverse_smith_waterman_sse2(ref, segLen, 2, 1, 2, 1, queryProfile, end_ref, &begin_ref, &begin_read, end_seg, score, 4);
	queryProfile_destructor(queryProfile);
	
	// For Banded SmithWaterman:
	pair< pair<unsigned int, unsigned int>, pair<unsigned int, unsigned int> > hr;
	hr.first.first   = 0;                     // reference begin
	hr.first.second  = end_ref - begin_ref;   // reference end
	hr.second.first  = 0;                     // query begin
	hr.second.second = end_read - begin_read; // query end
	
	const unsigned int referenceLen = strlen(ref) - begin_ref + 1;
	const unsigned int queryLen     = strlen(read) - begin_read + 1;
	const unsigned char bandwidth = abs(int(referenceLen - queryLen)) + 7;
	CBandedSmithWaterman bsw(2, -1, 2, 1, bandwidth);
	const char* pReference = &ref[begin_ref - 1];
	const char* pQuery = &read[begin_read - 1];
	unsigned int referenceBSW;
	string cigarBSW;
	bsw.Align(referenceBSW, cigarBSW, pReference, referenceLen, pQuery, queryLen, hr);
	printf("reference:    %s %3u\n", cigarBSW.c_str(), referenceBSW);	*/
	fprintf(stdout, "max score: %d, end_ref: %d, end_read: %d\nmax2 score: %d, end_ref: %d, end_read: %d\n", 
			bests[0].score, bests[0].ref, bests[0].read, bests[1].score, bests[1].ref, bests[1].read);
	free(ref);
	free(read);
    return 0;
}



/*
 *  ssw.cpp
 *
 *  Created by Mengyao Zhao on 04/05/11.
 *	Version 0.1.3
 *	New features: Record the highest 2*read_length scores among the highest score of each reference position. 
 *	Give out the most possible 2nd distinguished best alignment score as well as the best alignment score and
 *	its ending position. Passed the test with real data.
 *
 */

using namespace std;

#include <stdlib.h>
#include <vector>
#include <emmintrin.h>
#include "ssw.h"
#include "Benchmark.h"
#define SIZE 64

//Read reference and read sequences from the test file.
char* readFasta (vector<char*>* reads, char* file_ref, char* file_read) {
	//Read reference sequence.
	char* ref;
	FILE* Handle = fopen (file_ref, "r");
	if(Handle == NULL) {
		fprintf(stderr, "Error: Can't access reference file: %s.\n", file_ref);
	}
	unsigned int s = SIZE, i = 0;
	char* line = (char*) calloc(SIZE, 1);
	while (! feof (Handle)) {
		char letter = fgetc (Handle);
		if (letter == '\n') {
			line[i] = '\0';
			if (!strstr(line, ">") && i != 0) {
				ref = (char*) calloc (i + 1, 1);
				strcpy (ref, line);
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
	
	//Read reads.
	Handle = fopen (file_read, "r");
	if(Handle == NULL) {
		fprintf(stderr, "Error: Can't access read file: %s.\n", file_read);
	}
	s = SIZE, i = 0;
	while (! feof (Handle)) {
		char letter = fgetc (Handle);
		if (letter == '\n') {
			line[i] = '\0';
			if (!strstr(line, ">") && i != 0) {
				char* temp = (char*) calloc(i + 1, 1);
				strcpy (temp, line);
				reads->push_back(temp);
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
	return ref;
}

int main (int argc, char * const argv[]) {
	// validate argument count
	if (argc < 2) {
		fprintf (stderr, "USAGE: %s <ref.fa> <reads.fa>\n", argv[0]);
		exit (1);
	}
	
	char* file_ref = argv[1];
	char* file_read = argv[2];
	vector<char*> reads;
	char* ref = readFasta(&reads, file_ref, file_read);
	unsigned int refLen = strlen(ref); 
	unsigned int* ref_num = ref_amino2num (ref, refLen);
	
	// start timing the algorithm
	CBenchmark bench;
	bench.Start();	
	for (vector<char*>::iterator readIt = reads.begin(); readIt < reads.end(); readIt ++) {
		fprintf(stderr, "%s\n", (*readIt));
		char* read = (char*) calloc(strlen(*readIt) + 1, 1);
		strcpy (read, (*readIt));
		
		__m128i* vProfile = queryProfile_constructor(read, 2, 1, 4);
		unsigned int end_seg = 0;
		unsigned int readLen = strlen(read);
		alignment_end* bests = smith_waterman_sse2(ref_num, refLen, readLen, 2, 1, 2, 1, vProfile, &end_seg, 4);
		free(vProfile);
		
		if (bests[0].score != 0) {
			fprintf(stdout, "max score: %d, end_ref: %d\nmax2 score: %d, end_ref: %d\n", 
					bests[0].score, bests[0].ref, bests[1].score, bests[1].ref);		
			}else {
			fprintf(stdout, "No alignment found for this read.\n");
		}
		free(read);
	}
	free(ref);
	
	// stop timing the algorithm
	bench.Stop();
	bench.DisplayTime("Striped SmithWatterman");
    return 0;
}


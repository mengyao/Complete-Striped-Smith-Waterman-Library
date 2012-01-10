/*  main.c
 *  Created by Mengyao Zhao on 06/23/11.
 *	Version 0.1.4
 *  Last revision by Mengyao Zhao on 01/10/12.
 *	New features: find the best alignment beginning position 
 */

#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <zlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "ssw.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

char* seq_reverse(const char* seq, int32_t end)	/* end is 0-based alignment ending position */	
{									
	char* reverse = (char*)calloc(end + 2, sizeof(char*));	
	int32_t start = 0;
	reverse[end + 1] = '\0';				
	while (start <= end) {			
		reverse[start] = seq[end];		
		reverse[end] = seq[start];		
		++ start;					
		-- end;						
	}								
	return reverse;					
}									

int main (int argc, char * const argv[]) {
	// validate argument count
	if (argc < 2) {
		fprintf (stderr, "USAGE: %s <ref.fa/fq> <reads.fa/fq>\n", argv[0]);
		exit (1);
	}

	clock_t start, end;
	float cpu_time;
	gzFile ref_fp;
	kseq_t *ref_seq;
	int l;
	int m;
	ref_fp = gzopen(argv[1], "r");
	ref_seq = kseq_init(ref_fp);

	start = clock();
	while ((l = kseq_read(ref_seq)) >= 0) {
		printf("ref_name: %s\n", ref_seq->name.s);
		gzFile read_fp = gzopen(argv[2], "r");
		kseq_t*	read_seq = kseq_init(read_fp);
		int32_t refLen = strlen(ref_seq->seq.s);
		char* ref_reverse = seq_reverse(ref_seq->seq.s, refLen - 1); 
		fprintf(stderr, "reverse_ref: %s\n", ref_reverse); 										
		while ((m = kseq_read(read_seq)) >= 0) {
			char *read_reverse;
			alignment_end *bests, *bests_reverse;
			printf("read_name: %s\n", read_seq->name.s);
			printf("read_seq: %s\n", read_seq->seq.s); 
			
			int32_t readLen = strlen(read_seq->seq.s);
			__m128i* vProfile = queryProfile_constructor(read_seq->seq.s, 2, 1, 4);
	//		__m128i* vProfile = queryProfile_constructor(read_seq->seq.s, 2, 2, 4);
			bests = smith_waterman_sse2(ref_seq->seq.s, refLen, readLen, 2, 1, 2, 1, vProfile, 0, 4);
	//		alignment_end* bests = smith_waterman_sse2(ref_seq->seq.s, refLen, readLen, 3, 1, 3, 1, vProfile, 0, 4);
			free(vProfile);
			
			read_reverse = seq_reverse(read_seq->seq.s, bests[0].read);
			//ref_reverse = seq_reverse(ref_seq->seq.s, bests[0].ref);
			fprintf(stderr, "reverse_read: %s\n", read_reverse); 										
			vProfile = queryProfile_constructor(read_reverse, 2, 1, 4);
			bests_reverse = smith_waterman_sse2(ref_reverse + refLen - bests[0].ref - 1, bests[0].ref + 1, bests[0].read + 1, 2, 1, 2, 1, vProfile, bests[0].score, 4);
			free(vProfile);
			free(read_reverse);
			
			if (bests[0].score != 0) {
				fprintf(stdout, "max score: %d, end_ref: %d, end_read: %d\nbegin_ref: %d, begin_read: %d\n", 
						bests[0].score, bests[0].ref + 1, bests[0].read + 1, bests[0].ref - bests_reverse[0].ref + 1, bests[0].read - bests_reverse[0].read + 1);		
				}else {
				fprintf(stdout, "No alignment found for this read.\n");
			}
		}
		free(ref_reverse);
		kseq_destroy(read_seq);
		gzclose(read_fp);
	}
	end = clock();
	cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
	fprintf(stdout, "CPU time: %f seconds\n", cpu_time);

	kseq_destroy(ref_seq);
	gzclose(ref_fp);
	return 0;
}
